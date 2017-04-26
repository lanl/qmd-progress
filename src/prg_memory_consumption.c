#include <errno.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>

#define BUFFER_LENGTH 1000

void prg_memory_consumption(long long int *vm_peak, long long int *vm_size,
        long long int *pid)
{
    FILE *status;
    char buffer[BUFFER_LENGTH];
    char *substring;

    if ((status = fopen("/proc/self/status", "r")) == NULL)
    {
        fprintf(stderr, "error reading memory consumption: %s\n", strerror(errno));
        exit(1);
    }

    *vm_peak = -1;
    *vm_size = -1;
    *pid = -1;

    while (fgets(buffer, BUFFER_LENGTH, status))
    {
        if (strstr(buffer, "Pid")) {
            strtok(buffer, " \t");
            *pid = strtoll(strtok(NULL, " "), NULL, 10);
        }

        if (strstr(buffer, "VmPeak")) {
            strtok(buffer, " \t");
            *vm_peak = strtoll(strtok(NULL, " "), NULL, 10) / 1024;
        }

        if (strstr(buffer, "VmSize")) {
            strtok(buffer, " \t");
            *vm_size = strtoll(strtok(NULL, " "), NULL, 10) / 1024;
        }
    }

    fclose(status);

    if (*vm_peak < 0 || *vm_size < 0)
    {
        fprintf(stderr, "could not read memory consumption\n");
        exit(1);
    }
}
