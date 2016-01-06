#include <iostream>
#include <fstream>
#include <mpi.h>
#include <vector>
#include <cmath>
#include <iomanip>
#include <sstream>

std::ostream &debug(int rank) {
#ifdef DEBUG
    return std::cerr << "[" << rank << "] ";
#else
    static std::stringstream ss;
    return ss;
#endif
}

class Field
{
    MPI_Comm comm;
    int rank, nproc, top, bottom;

    int n, m;
    double *data, *prev_data, *buff;
    bool is_data_inited;

    double width, height;
    double dx, dy, dt, dxtam2, dytam2;
    double time;

    void init_field() {
        // Размер сетки выбираем исходя из заданных размеров задачи и точности.
        int all_n = (int)std::ceil(width / dx),
            all_m = (int)std::ceil(height / dy / nproc) * nproc;

        // И уточняем dx и dy. Они не будут больше, чем указано в точности. Т.е. точность может только повысится.
        dx = (double)width / all_n;
        dy = (double)height / all_m;

        // На каждом процессе будет вычисляться только all_m / nproc строк.
        n = all_n;
        m = all_m / nproc;

        debug(rank) << "Field init: " << m << " x " << n << " of " << all_m << " x " << all_n << "\n";

        // Создаем массив длинны width x height + 2 строки (соседей).
        const int datasize = n * (m + 2);
        data = new double[datasize];
        prev_data = new double[datasize];
        buff = new double[datasize];

        // Подобрал коэффициент alpha. Просчитаем также множители, чтобы уменьшить число операций в численной схеме.
        const double alpha = 1.0 / 4200.0;
        dxtam2 = dt * alpha / (dx * dx);
        dytam2 = dt * alpha / (dy * dy);
    }

    void init_data() {
        // Мы только что заварили чай в очень длинной паралелепипедо-образной таре
        // и рассматриваем температуру в срезе посередине длины.
        const double init_temp = 100.0;

        for (int i = 0; i < m + 2; ++i) {
            for (int j = 0; j < n; ++j) {
                data[i * n + j] = init_temp;
            }
        }
        is_data_inited = true;
    }

    void update_borders() {
        // Рассматриать будем, что температура среды не меняется.
        // Сверху мы будем дуть на чай, температура там пусть будет 20, температура комнаты.
        // Стенки справа и слева пусть будут температуры повыше, 40.
        // Стенка снизу, что стоит на полу будут ещё выше, 50.
        // Такая задача очень искусственная и захардкожена.

        for (int i = 0; i < m + 2; ++i) {
            data[i * n + 0] = data[i * n + (n - 1)] = 40;
        }

        if (top == MPI_PROC_NULL) {
            for (int j = 0; j < n; ++j) {
                data[0 * n + j] = 20;
            }
        }
        if (bottom == MPI_PROC_NULL) {
            for (int j = 0; j < n; ++j) {
                data[(m + 1) * n + j] = 50;
            }
        }
    }

    void calc_layer() {
        std::swap(data, prev_data);

        // Явная разностная схема, шаблон "крест".
        // y(i,j,t+1) = y(i,j,t) + dt / 4200.0 * (
        //     (y(i-1,j,t) + y(i+1,j,t) - 2 * y(i,j,t)) / (dx * dx) +
        //     (y(i,j-1,t) + y(i,j+1,t) - 2 * y(i,j,t)) / (dy * dy) )
        for (int i = 1; i < m + 1; ++i) {
            for (int j = 1; j < n - 1; ++j) {
                data[i * n + j] = prev_data[i * n + j] + (
                    dxtam2 * (prev_data[i * n + j - 1] + prev_data[i * n + j + 1] - 2 * prev_data[i * n + j]) +
                    dytam2 * (prev_data[(i - 1) * n + j] + prev_data[(i + 1) * n + j] - 2 * prev_data[i * n + j])
                );
            }
        }
    }

    // Пересылка краевых значений между соседями.
    void sendrecv_borders() {
        MPI_Sendrecv(data + m * n, n, MPI_DOUBLE, bottom, 0,
                     data, n, MPI_DOUBLE, top, 0, comm, MPI_STATUS_IGNORE);
        MPI_Sendrecv(data + n, n, MPI_DOUBLE, top, 0,
                     data + (m + 1) * n, n, MPI_DOUBLE, bottom, 0, comm, MPI_STATUS_IGNORE);
    }

public:
    Field(double h, double w, double _dx, double _dy, double _dt)
        : width(w), height(h), dx(_dx), dy(_dy), dt(_dt)
    {
        // Создаем коммуникатор с 1-мерной топологией сетки без цикличности.
        MPI_Comm_size(MPI_COMM_WORLD, &nproc);
        int periods[] = { 0 };
        MPI_Cart_create(MPI_COMM_WORLD, 1, &nproc, periods, 1, &comm);
        // Узнаем свой ранг и ранг соседей.
        MPI_Comm_rank(comm, &rank);
        MPI_Cart_shift(comm, 0, 1, &top, &bottom);

        init_field();
    }

    ~Field() {
        MPI_Comm_free(&comm);
        delete[] data;
        delete[] prev_data;
        delete[] buff;
    }

    // Просчитывает следующие _time секунд.
    void calc_time(double _time, bool fromstart = false) {
        if (!is_data_inited || fromstart) {
            init_data();
            time = 0;
        }
        for (double t = 0; t < _time; t += dt, time += dt) {
            update_borders();
            calc_layer();
            sendrecv_borders();
        }
    }

    // Вся матрица дописывается в конец файла.
    void print_to_file(const char *filename) {
        if (rank == 0) {
            std::ofstream fout(filename, std::ios::app);

            fout << time << "\n" << m * nproc << " " << n << "\n";
            for (int i = 1; i < m + 1; ++i) {
                for (int j = 0; j < n; ++j) {
                    fout << std::fixed << std::setw(11) << std::setprecision(6)  << data[i * n + j] << " ";
                }
                fout << "\n";
            }
            for (int k = 1; k < nproc; ++k) {
                MPI_Recv(buff, m * n, MPI_DOUBLE, k, 1, comm, MPI_STATUS_IGNORE);
                for (int i = 0; i < m; ++i) {
                    for (int j = 0; j < n; ++j) {
                        fout << std::fixed << std::setw(11) << std::setprecision(6) << buff[i * n + j] << " ";
                    }
                    fout << "\n";
                }
            }

            fout.close();
        } else {
            MPI_Send(data + n, m * n, MPI_DOUBLE, 0, 1, comm);
        }
    }
};

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

    wait_for_debug();

    {
        int rank;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);

        const char *fname = "info.txt";

        // Затираем файл.
        std::ofstream ff(fname);
        ff.close();

        // Много хардкода, не хотелось тут заморачиваться с конфигом.
        Field field(1, 1, 0.01, 0.01, 0.01);

        // Записываем в файл первые 300 секунд симуляции.
        for (int i = 0; i < 300; ++i) {
            field.calc_time(1);
            field.print_to_file(fname);

            if (rank == 0) {
                debug(rank) << i << "/" << 200 << "\n";
            }
        }
    }

    MPI_Finalize();
    return 0;
}
