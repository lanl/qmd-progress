#include <errno.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>

#define BUFFER_LENGTH 1000

void
prg_memory_consumption(
    long long int *vm_peak,
    long long int *vm_size,
    long long int *pid,
    long long int *ppid)
{
    FILE *status;
    char buffer[BUFFER_LENGTH];
    char *substring;

    if ((status = fopen("/proc/self/status", "r")) == NULL)
    {
        fprintf(stderr, "error reading memory consumption: %s\n",
                strerror(errno));
        exit(1);
    }

    *vm_peak = -1;
    *vm_size = -1;
    *pid = -1;
    *ppid = -1;

    while (fgets(buffer, BUFFER_LENGTH, status))
    {
        substring = strtok(buffer, " \t");
        if (!substring)
            continue;
        if (strcmp(substring, "Pid:") == 0)
        {
            *pid = strtoll(strtok(NULL, " \t"), NULL, 10);
        }

        if (strcmp(substring, "PPid:") == 0)
        {
            *ppid = strtoll(strtok(NULL, " \t"), NULL, 10);
        }

        if (strcmp(substring, "VmPeak:") == 0)
        {
            *vm_peak = strtoll(strtok(NULL, " \t"), NULL, 10) / 1024;
        }

        if (strcmp(substring, "VmSize:") == 0)
        {
            *vm_size = strtoll(strtok(NULL, " \t"), NULL, 10) / 1024;
        }
    }

    fclose(status);

    if (*vm_peak < 0 || *vm_size < 0)
    {
        fprintf(stderr, "could not read memory consumption\n");
        exit(1);
    }
}
