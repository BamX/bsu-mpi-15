#include <iostream>
#include <mpi.h>
#include <fstream>
#include <cmath>
#include <sstream>

std::ostream &debug(int rank) {
#ifdef DEBUG
    return std::cerr << "[" << rank << "] ";
#else
    static std::stringstream ss;
    return ss;
#endif
}

class MatrixMult
{
    int q, n, coords[2], root, cart_root;
    double *pA, *pcA, *pB, *pC;

    int top, bottom, rank, cart_rank;
    MPI_Comm comm, commRow, commCol;
    MPI_Datatype blockType;

    // 1-й этап. Все процессы в строке получают блок A из столбца root_col.
    void bcast_block_a(int root_col) {
        if (coords[1] == root_col) {
            memcpy(pcA, pA, n * n * sizeof(double));
        }
        MPI_Bcast(pcA, n * n, MPI_DOUBLE, root_col, commRow);
    }

    // 2-й этап. Перемножение блоков.
    void multiply_blocks() {
        for (int i = 0; i < n; ++i) {
            for (int k = 0; k < n; ++k) {
                for (int j = 0; j < n; ++j) {
                    pC[i * n + j] += pcA[i * n + k] * pB[k * n + j];
                }
            }
        }
    }

    // 3-й этап. Циклическая отправка блока B вверх по столбцу.
    void sendrecv_block_b(int root_col) {
        MPI_Sendrecv_replace(pB, n * n, MPI_DOUBLE, top, root_col, bottom, root_col, commCol, MPI_STATUS_IGNORE);
    }

    // Определяем процесс, который будет считывать и записывать матрицу
    // Этот процесс должен иметь доступ к файлу.
    // Из всех таких процессов выбирается с наименьшим рангом.
    void find_root_by_filename(const char *filename, MPI_Comm scomm) {
        std::ifstream fcheck(filename);
        int fileExists = fcheck.good() ? 1 : 0;
        fcheck.close();

        int size, rank;
        MPI_Comm_size(scomm, &size);
        MPI_Comm_rank(scomm, &rank);

        int *processes = new int[size];
        MPI_Allgather(&fileExists, 1, MPI_INT, processes, 1, MPI_INT, scomm);

        root = MPI_PROC_NULL;
        for (int i = 0; i < size; ++i) {
            if (processes[i] == 1) {
                root = i;
            }
        }

        delete[] processes;

        if (root == MPI_PROC_NULL) {
            if (rank == 0) {
                std::cerr << "File not found from any process: " << filename << "\n";
            }
            MPI_Finalize();
            exit(1);
        }
    }

    // root считывает всю матрицу из fin и раздает всем процессам блоки в массив part.
    void read_and_scatter_matrix(std::ifstream &fin, double *part) {
        MPI_Request recvRequest;
        MPI_Irecv(part, n * n, MPI_DOUBLE, cart_root, 0, comm, &recvRequest);

        if (rank == root) {
            double *matrix = new double[(n * q) * n];

            for (int k = 0; k < q; ++k) {
                int idx = 0;
                for (int i = 0; i < n; ++i) {
                    for (int j = 0; j < n * q; ++j) {
                        fin >> matrix[idx++];
                    }
                }

                for (int j = 0; j < q; ++j) {
                    int rank, coords[] = { k, j };
                    MPI_Cart_rank(comm, coords, &rank);
                    MPI_Send(matrix + j * n, 1, blockType, rank, 0, comm);
                }
            }

            delete[] matrix;
        }

        MPI_Wait(&recvRequest, MPI_STATUS_IGNORE);
    }

    // Со всех процессов в root собираются блоки матрицы C и записываются в файл.
    void gather_and_save_matrix(const char *filename) {
        double *matrix = NULL;

        MPI_Request sendRequest;
        MPI_Isend(pC, n * n, MPI_DOUBLE, cart_root, coords[0], comm, &sendRequest);

        if (rank == root) {
            matrix = new double[(n * q) * n];

            std::ofstream fout(filename);
            fout << n * q << "\n";

            for (int i = 0; i < q; ++i) {
                // Для каждой строки процессов собираем блоки
                for (int j = 0; j < q; ++j) {
                    MPI_Status st;
                    MPI_Probe(MPI_ANY_SOURCE, i, comm, &st);
                    int rank = st.MPI_SOURCE, coords[2];
                    MPI_Cart_coords(comm, rank, 2, coords);
                    MPI_Recv(matrix + coords[1] * n, 1, blockType, rank, i, comm, MPI_STATUS_IGNORE);
                }

                int idx = 0;
                for (int i = 0; i < n; ++i) {
                    fout << matrix[idx++];
                    for (int j = 1; j < n * q; ++j) {
                        fout << " " << matrix[idx++];
                    }
                    fout << "\n";
                }
            }

            fout.close();
            delete[] matrix;
        }

        MPI_Wait(&sendRequest, MPI_STATUS_IGNORE);
    }

    // Отладочный вывод блока ar.
    void print(double *ar) {
#ifdef DEBUG
        for (int ii = 0; ii < q; ++ii) {
            for (int jj = 0; jj < q; ++jj) {
                if (coords[0] == ii && coords[1] == jj) {
                    std::cerr << coords[0] << " " << coords[1] << "\n";
                    for (int i = 0; i < n; ++i) {
                        for (int j = 0; j < n; ++j) {
                            std::cerr << ar[i * n + j] << " ";
                        }
                        std::cerr << "\n";
                    }
                }
                MPI_Barrier(comm);
            }
        }
#endif
    }

public:
    MatrixMult(const char *filename, const MPI_Comm &_comm)
    {
        MPI_Comm scomm;
        MPI_Comm_dup(_comm, &scomm);

        find_root_by_filename(filename, scomm);

        int size;
        MPI_Comm_size(scomm, &size);
        MPI_Comm_rank(scomm, &rank);

        // Проверяем, что число процессов является квадратом.
        q = (int)std::sqrt(size);
        if (q * q != size) {
            std::cerr << "Processes count have to be square\n";
            MPI_Finalize();
            exit(1);
        }

        int all_n;
        std::ifstream fin;

        if (rank == root) {
            fin = std::ifstream(filename);
            fin >> all_n;
        }

        MPI_Bcast(&all_n, 1, MPI_INT, root, scomm);

        // Размер матрици должен быть кратен корню числа процессов.
        if (all_n % q != 0) {
            std::cerr << "Matrix wrong size\n";
            MPI_Finalize();
            exit(2);
        }

        // Вычисляем размер блока и создаем коммуникатор с топологией 2-мерной сетки q x q.
        n = all_n / q;
        int dims[] = { q, q };
        int periods[] = { 1, 1 };
        MPI_Cart_create(scomm, 2, dims, periods, 1, &comm);

        // Также создаем коммуникатор для строки в сетке.
        int remain_row[] = { 0, 1 }, remain_col[] = { 1, 0 };
        MPI_Cart_sub(comm, remain_row, &commRow);
        MPI_Cart_sub(comm, remain_col, &commCol);

        // Получаем наши координаты в сетке.
        MPI_Comm_rank(comm, &cart_rank);
        MPI_Cart_coords(comm, cart_rank, 2, coords);

        // Синхронизируем root в сетке.
        cart_root = cart_rank;
        MPI_Bcast(&cart_root, 1, MPI_INT, root, scomm);

        // Получаем соседей выше и ниже.
        MPI_Cart_shift(commCol, 0, 1, &top, &bottom);

        // Создаем тип для блока матрицы.
        MPI_Type_vector(n, n, n * q, MPI_DOUBLE, &blockType);
        MPI_Type_commit(&blockType);

        pA = new double[n * n];
        pcA = new double[n * n];
        pB = new double[n * n];
        pC = new double[n * n];

        read_and_scatter_matrix(fin, pA);
        read_and_scatter_matrix(fin, pB);
        memset(pC, 0, n * n * sizeof(double));

        if (fin.is_open()) {
            fin.close();
        }

        MPI_Comm_free(&scomm);
    }

    ~MatrixMult() {
        MPI_Comm_free(&comm);
        MPI_Comm_free(&commRow);
        delete[] pA;
        delete[] pcA;
        delete[] pB;
        delete[] pC;
    }

    void multiply_and_save(const char *filename) {
        for (int i = 0; i < q; ++i) {
            bcast_block_a((coords[0] + i) % q);
            multiply_blocks();
            sendrecv_block_b(i);
        }

        gather_and_save_matrix(filename);
    }
};

void generate_matrix(const char *filename, int n) {
    std::ofstream fout(filename);

    fout << n << "\n";

    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            fout << (i == j ? 1 : 0) << " ";
        }
        fout << "\n";
    }

    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            fout << (i == j ? 2 : 0) << " ";
        }
        fout << "\n";
    }

    fout.close();
}

// Функция для отладки в IDE. Дает возможность подключиться к процессу.
void wait_for_debug() {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if (rank == 0) {
        int i = 0;
        while (i == 0) {
            sleep(0);
        }
    }
    MPI_Barrier(MPI_COMM_WORLD);
}

int main(int argc, char **argv) {
    MPI_Init(&argc, &argv);

    //wait_for_debug();

    const char *ifname = "matrix.txt";
    const char *ofname = "out.txt";
    generate_matrix(ifname, 12);

    {
        MatrixMult mm(ifname, MPI_COMM_WORLD);
        mm.multiply_and_save(ofname);
    }
    
    MPI_Finalize();
    return 0;
}
