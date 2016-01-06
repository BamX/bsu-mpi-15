#include <iostream>
#include <sstream>
#include <string>

#ifndef TASK
#define TASK 1
#endif

std::ostream &debug() {
#ifdef DEBUG
    return std::cerr;
#else
    static std::stringstream ss;
    return ss;
#endif
}

namespace shared {

    // Ограничение максимального числа процессов
    static int const MAX_P = 30;
    // Шаг интегрирования
    static double const INTG_STEP = 1e-8;

    // Число процессов
    int n;
    // Границы и длина части рассчета
    double lx, rx, dx;

    // Целевая ф-я
    double foo(double x) {
        return x * (x + 1);
    }

    // Сообщение об ошибке
    void usage() {
        std::cerr << "Usage: app [proc_count] [lx] [rx]\n";
        exit(1);
    }

    // Получение параметров из аргументов запуска
    void parseArgs(int argc, const char * argv[]) {
        if (argc != 4) usage();

        n = atoi(argv[1]);
        if (n <= 0 || n >= MAX_P) usage();

        lx = atof(argv[2]);
        rx = atof(argv[3]);
        if (lx >= rx) usage();
        
        dx = (rx - lx) / n;
    }
}

#if TASK == 1

#include <semaphore.h>

namespace task01 {

    int rank = 0;
    int children[shared::MAX_P], pipes[2];
    double l, r;

    void makeForks() {
        l = shared::lx; r = l + shared::dx;

        // Открытие пайпов для сбора результатов
        pipe(pipes);
        for (int i = 0; i < shared::n - 1; ++i) {
            int chId = fork();

            // Один процесс порождает множество других.
            if (chId > 0) {  // parent
                children[i] = chId;
            } else { // child
                rank = i + 1;
                break;
            }

            l += shared::dx; r += shared::dx;
        }

        // 0-й процесс будет собирать данные
        if (rank == 0) {
            close(pipes[1]);
        } else {
            close(pipes[0]);
        }
    }

    // Вычисление интеграла на отрезке
    double calc() {
        double sum = 0;
        for (double t = l; t < r; t += shared::INTG_STEP) {
            sum += shared::foo(t) * shared::INTG_STEP;
        }
        return sum;
    }

    // Сбор данных
    double gather(double value) {
        // Семафор для синхронизации
        sem_t *sharedSem = sem_open("BSU_TASK01", O_CREAT, 0644, 1);

        if (rank == 0) {
            for (int i = 0; i < shared::n - 1; ++i) {
                double otherSum;
                while (read(pipes[0], &otherSum, sizeof(double) * 1) < 0) {
                    switch (errno) {
                        case EINTR:
                            std::cerr << "EINTR\n";
                            break;
                        case EAGAIN:
                            std::cerr << "EAGAIN\n";
                            break;
                        case EIO:
                            std::cerr << "EIO\n";
                            break;
                        case EISDIR:
                            std::cerr << "EISDIR\n";
                            break;
                        case EBADF:
                            std::cerr << "EBADF\n";
                            break;
                        case EINVAL:
                            std::cerr << "EINVAL\n";
                            break;
                        case EFAULT:
                            std::cerr << "EFAULT\n";
                            break;
                        default:
                            std::cerr << "UNKNOWN\n";
                            break;
                    }
                    sleep(0);
                }
                value += otherSum;
            }
        } else {
            sem_wait(sharedSem);
            write(pipes[1], &value, sizeof(double) * 1);
            sem_post(sharedSem);
        }

        sem_close(sharedSem);
        sem_unlink("BSU_TASK01");
        return value;
    }

    // Завершение дочерних процессов
    void join() {
        if (rank == 0) {
            for (int i = 0; i < shared::n - 1; ++i) {
                waitpid(children[i], NULL, 0);
            }
        }
    }

    // Основная функция решения задачи
    void process() {
        makeForks();
        double value = calc();
        value = gather(value);
        join();

        if (rank == 0) {
            printf("%.6f\n", value);
        }
    }
}

#elif TASK == 2

#include <thread>
#include <mutex>
#include <queue>

namespace task02 {

    // Потоки обработки
    std::thread *threads[shared::MAX_P];

    // Мьютекс синхронизации очереди обработки отрезков вычисления
    std::mutex queueMutex;
    std::queue<double> queue;

    // Мьютекс формирования результата
    std::mutex valueMutex;
    double value;

    // Получение левой границы отрезка из очереди
    bool getLeft(double *l) {
        std::lock_guard<std::mutex> lock(queueMutex);
        if (queue.empty()) return false;
        *l = queue.front();
        queue.pop();
        return true;
    }

    // Формирование результата
    void addValue(double x) {
        std::lock_guard<std::mutex> lock(valueMutex);
        value += x;
    }

    // Функция, выполняемая в каждом потоке для просчета результата
    void calcSync() {
        double l;
        if (!getLeft(&l)) return;

        double sum = 0;
        for (double t = l, r = l + shared::dx; t < r; t += shared::INTG_STEP) {
            sum += shared::foo(t) * shared::INTG_STEP;
        }
        addValue(sum);
    }

    // Общая функция вычисления интеграла на всём отрезке
    void calcAsync() {
        value = 0;

        // Подотрезки добавляются в очередь обработки
        for (int i = 0; i < shared::n; ++i) {
            queue.push(shared::lx + i * shared::dx);
        }

        // Запускаются потоки обработки (на 1 меньше n, т.к. есть и главный поток)
        for (int i = 0; i < shared::n - 1; ++i) {
            threads[i] = new std::thread(calcSync);
        }

        // Происходит вычисление на главном потоке
        calcSync();

        // Джоин потоков.
        for (int i = 0; i < shared::n - 1; ++i) {
            threads[i]->join();
            delete threads[i];
        }
    }
    
    void process() {
        calcAsync();
        printf("%.6f\n", value);
    }
}

#endif

int main(int argc, const char * argv[]) {
    shared::parseArgs(argc, argv);
#if TASK == 1
    task01::process();
#elif TASK == 2
    task02::process();
#endif
    return 0;
}
