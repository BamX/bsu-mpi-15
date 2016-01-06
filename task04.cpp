#include <iostream>
#include <mpi.h>

// Создает коммуникатор с топологией "бинарное дерево" с корнем в процессе root.
// Создание коммуникатора происходит при помощи функции MPI_Graph_create.
// Полученный граф(дерево) является ориентированным.
int MPI_Binary_tree_create(MPI_Comm comm_old, int root, int reorder, MPI_Comm *comm) {
    int size;
    MPI_Comm_size(comm_old, &size);

    int *index = new int[size], *edges = new int[size * 2];
    int edges_count = 0;
    // Для каждого процесса в коммуникаторе определяем до 2-х потомков.
    for (int i = 0; i < size; ++i) {
        // Определяем номер вершины в дереве(вершины нумеруются сверху вниз слева направо).
        int node_id = (i - root + size) % size;
        for (int j = 0; j < 2; ++j) {
            // Определяем номер вершины потомка и проверяем, попадает ли номер потомка в число процессов.
            int chld_shift = node_id * 2 + (j + 1);
            if (chld_shift < size) {
                edges[edges_count++] = (root + chld_shift) % size;
            }
        }
        index[i] = edges_count;
    }

    int retval = MPI_Graph_create(comm_old, size, index, edges, reorder, comm);
    delete[] index;
    delete[] edges;
    return retval;
}

// Аналог MPI_Bcast, для обмена данными создается новый коммуникатор с топологией "бинарное дерево".
int MPI_Fbcast(void *buffer, int count, MPI_Datatype datatype, int root, MPI_Comm comm) {
    MPI_Comm tree_comm;
    MPI_Request requests[2];
    int rank, tree_rank;
    int children[2], children_count;

    // Создаем коммуникатор требуемой топологии.
    MPI_Binary_tree_create(comm, root, 1, &tree_comm);

    // Если мы не являемся корнем bcast'а, ждем и получаем данные от предка в дереве.
    MPI_Comm_rank(comm, &rank);
    if (rank != root) {
        MPI_Recv(buffer, count, datatype, MPI_ANY_SOURCE, 0, tree_comm, MPI_STATUS_IGNORE);
    }

    // Получаем информацию о потомках в дереве.
    MPI_Comm_rank(tree_comm, &tree_rank);
    MPI_Graph_neighbors_count(tree_comm, tree_rank, &children_count);

    // Если потомки есть, передаем им данные.
    if (children_count > 0) {
        MPI_Graph_neighbors(tree_comm, tree_rank, 2, children);

        for (int chid = 0; chid < children_count; ++chid) {
            MPI_Isend(buffer, count, datatype, children[chid], 0, tree_comm, &requests[chid]);
        }

        // Передача происходит асинхронно обоим потомкам, но перед выходом из bcast'a дожидаемся её завершения.
        MPI_Waitall(children_count, requests, MPI_STATUSES_IGNORE);
    }

    // Перед выходом также освобождаем временный коммуникатор.
    MPI_Comm_free(&tree_comm);
    return 0;
}

typedef int (*Bcast_func)(void *, int, MPI_Datatype, int, MPI_Comm);

const int root = 13;
const int repeats = 1000;
const int buf_size = 1 << 20;
char buf[buf_size];

size_t mesure_bcast(Bcast_func bcast, int repeats, int buf_size) {
    auto start_time = std::chrono::high_resolution_clock::now();

    for (int i = 0; i < repeats; ++i)
        bcast(buf, buf_size, MPI_CHAR, root, MPI_COMM_WORLD);

    auto end_time = std::chrono::high_resolution_clock::now();

    return std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time).count() / repeats;
}

int main(int argc, char **argv) {
    MPI_Init(&argc, &argv);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    if (rank == root) {
        sprintf(buf, "Hello, world!");
    }

    if (rank == root) {
        std::cerr << "bytes\tMPI_Fbcast\tMPI_Bcast\n";
    }

    for (int i = 10; i < 21; ++i) {
        auto time_fbcast = mesure_bcast(MPI_Fbcast, repeats, 1 << i);
        auto time_bcast = mesure_bcast(MPI_Bcast, repeats, 1 << i);

        if (rank == root) {
            std::cerr << (1 << i) << "\t" << time_fbcast << "\t" << time_bcast << "\n";
        }
    }

    MPI_Finalize();
    return 0;
}
