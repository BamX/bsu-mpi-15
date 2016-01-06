#include <iostream>
#include <sstream>

#include <thread>
#include <mutex>
#include <semaphore.h>
#include <sys/ipc.h>
#include <sys/shm.h>

namespace task03 {

    static key_t const CONTROL_KEY = 13223; // ключ shared memory контрол-блока
    static int const ACCESS_MASK = 0666;
    static size_t const INITIAL_SIZE = 32;

    struct Control {
        key_t memKey; // текущий ключ shared memory блока сообщений
        size_t memSize; // capacity блока сообщений
        size_t currentSize; // фактический размер блока сообщений
        size_t refCount; // число активных доступов к контрол-блоку
    };

    // Именной семафор с сигнатурой для использования lock guard'а.
    class Semaphore {
        sem_t *sharedSem;
        const char *nameSem;
        bool haveToRemove;
        
    public:
        Semaphore (const char *name, int accessMask = 0644) : nameSem(name) {
            haveToRemove = false;
            sharedSem = sem_open(name, 0, accessMask, 1);
            if (sharedSem == SEM_FAILED) {
                sharedSem = sem_open(name, O_CREAT, accessMask, 1);
                sem_post(sharedSem);
                haveToRemove = true;
            }
        }

        ~Semaphore () {
            if (haveToRemove) {
                sem_unlink(nameSem);
            }
            sem_close(sharedSem);
        }

        void lock() {
            sem_wait(sharedSem);
        }

        void unlock() {
            sem_post(sharedSem);
        }
    };

    // Класс, реализующий список сообщений в shared memory
    class SharedMessages {
        // Семафор синхронизации изменения данных
        Semaphore sem;

        key_t localKey;
        int controlMem, messagesMem;

        // Область в shared memory с контрольной информацией
        Control *control;
        // Область в shared memory со списком сообщений разделенных \0
        char *buffer;

        // Инициализация колнтольного объекта
        void initControl() {
            controlMem = shmget(CONTROL_KEY, sizeof(Control), ACCESS_MASK);
            if (controlMem == -1) {
                controlMem = shmget(CONTROL_KEY, sizeof(Control), ACCESS_MASK | IPC_CREAT);
                control = (Control *)shmat(controlMem, NULL, 0);
                control->memKey = localKey = -1;
                control->memSize = control->currentSize = 0;
                control->refCount = 1;
            } else {
                control = (Control *)shmat(controlMem, NULL, 0);
                ++(control->refCount);
            }
        }

        // Инициализация списка сообщений
        void initMem() {
            if (control->memSize == 0) {
                enlargeMem();
            } else {
                checkMem();
            }
        }

        // Проверка, не изменился ли блок с сообщениями в shared memory.
        void checkMem() {
            if (control->memKey != localKey) {
                localKey = control->memKey;
                if (buffer != NULL) {
                    shmdt(buffer);
                }
                messagesMem = shmget(localKey, control->memSize, ACCESS_MASK);
                buffer = (char *)shmat(messagesMem, NULL, 0);
            }
        }

        // Расширение блока памяти сообщений при истечении текущего.
        // Расширение в 2 раза.
        void enlargeMem() {
            size_t newSize = std::max(INITIAL_SIZE * sizeof(char), control->memSize * 2);
            key_t newKey = (control->memKey == -1 ? CONTROL_KEY : control->memKey) + 1;
            int newMem = shmget(newKey, newSize, 0666 | IPC_CREAT);
            char *newBuf = (char *)shmat(newMem, NULL, 0);

            if (control->memSize > 0) {
                memcpy(newBuf, buffer, control->memSize);
                shmdt(buffer);
                shmctl(messagesMem, IPC_RMID, NULL);
            }

            buffer = newBuf;
            messagesMem = newMem;
            localKey = control->memKey = newKey;
            control->memSize = newSize;
        }

        // Получение размера блока сообщений при смещении fromByte байт.
        size_t getSize(size_t fromByte = 0) {
            if (control->currentSize < fromByte)
                return 0;

            return (control->currentSize - fromByte) / sizeof(char);
        }

    public:
        SharedMessages() : sem("BSU_TASK03") {
            std::lock_guard<Semaphore> lock(sem);

            localKey = -1;
            initControl();
            initMem();
        }
        ~SharedMessages() {
            std::lock_guard<Semaphore> lock(sem);

            checkMem();
            size_t rc = --control->refCount;
            shmdt(buffer);
            shmdt(control);
            if (rc == 0) {
                shmctl(controlMem, IPC_RMID, NULL);
                shmctl(messagesMem, IPC_RMID, NULL);
            }
        }

        // Дописывает сообщение в конец блока сообщений.
        // Если соробщение не помещается, происходит расширение.
        void pushMessage(const char *msg) {
            std::lock_guard<Semaphore> lock(sem);
            checkMem();

            size_t msgSize = strlen(msg) + 1;
            size_t newSize = control->currentSize + msgSize;
            while (control->memSize < newSize) {
                enlargeMem();
            }
            memcpy(buffer + control->currentSize, msg, msgSize);
            control->currentSize += msgSize;
        }

        // Получение блока сообщений при смещении fromByte.
        int getMessages(size_t *size, char **message, size_t fromByte = 0) {
            std::lock_guard<Semaphore> lock(sem);
            checkMem();

            *size = getSize(fromByte);
            if (*size > 0) {
                *message = new char[*size];
                memcpy(*message, buffer + fromByte, *size * sizeof(char));
                return 0;
            }
            return 1;
        }
    };

    SharedMessages messages;
    bool stop = false;

    // Функция запускается в отдельном потоке и выводит сообщения из messages.
    void printMessages() {
        size_t skipSize = 0;
        while (!stop) {
            size_t size;
            char *msg;

            if (messages.getMessages(&size, &msg, skipSize) == 0) {
                skipSize += size;

                char *output = msg;
                while (size > 0) {
                    std::cout << output << "\n";
                    size_t msgSize = strlen(output) + 1;
                    output += msgSize;
                    size -= msgSize;
                }

                delete[] msg;
            }

            std::this_thread::sleep_for(std::chrono::milliseconds(500));
        }
    }

    // Сервер просто выводит сообщения.
    void server() {
        std::cout << "Welcome to Server.\nPress any key to exit.\nMessages:\n";

        stop = false;
        std::thread printer(printMessages);
        getchar();
        stop = true;
        printer.join();
    }

    // Клиент принимает строки с сообщениями. Пустое сообщение приводит к выходу из системы.
    void client() {
        std::cout << "Welcome to Client.\nSend empty message to exit.\nYour messages:\n";

        while (true) {
            char buf[256];
            std::cin.getline(buf, 256, '\n');

            if (strlen(buf) == 0) break;

            messages.pushMessage(buf);
        }
    }

    // Вывод описания использования приложения.
    void usage() {
        std::cerr << "Usage:\n\tapp -c (for client)\n\tapp -s (for server)\n";
        exit(1);
    }

    // Запуск сервера/клиента.
    void process(int argc, const char * argv[]) {
        if (argc != 2) usage();

        if (strcmp(argv[1], "-s") == 0) {
            server();
        } else if (strcmp(argv[1], "-c") == 0) {
            client();
        } else {
            usage();
        }
    }
}

int main(int argc, const char * argv[]) {
    task03::process(argc, argv);
    return 0;
}
