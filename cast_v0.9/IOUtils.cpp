#ifndef IOUTILS_CPP
#define IOUTILS_CPP

#include <fstream>
#include <cstdlib>
#include <cstring>
#include <cerrno>

using namespace std;

void mustOpenStream(fstream &f, const char *filename, ios_base::openmode mode) {
  f.open(filename, mode);
  if (!f.is_open()) {
    fprintf(stderr, "Error -- failed to open file: %s\n", filename);
    exit(1);
  }
}

// from blatSrc/lib/common.c: mustOpen
FILE *mustOpenFile(const char *fileName, const char *mode)
/* Open a file - or squawk and die. */
{
FILE *f;

if (strcmp(fileName, "stdin") == 0)
    return stdin;
if (strcmp(fileName, "stdout") == 0)
    return stdout;
if ((f = fopen(fileName, mode)) == NULL)
    {
    //char *modeName = "";
    string modeName;
    if (mode)
        {
        if (mode[0] == 'r')
            modeName = " to read";
        else if (mode[0] == 'w')
            modeName = " to write";
        else if (mode[0] == 'a')
            modeName = " to append";
        }
    //errAbort("Can't open %s%s: %s", fileName, modeName, strerror(errno));
    fprintf(stderr, "Can't open %s%s: %s\n", fileName, modeName.c_str(), strerror(errno));
    }
return f;
}


#endif
