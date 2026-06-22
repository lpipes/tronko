/* This set of functions is designed for efficiently reading a large file
 * (either plaintext or .gz) line by line. While it is possible to do this by
 * repeatedly calling `gzgets`, this introduces the overhead of decompressing
 * every single line individually, which is costly when the lines are short. The
 * functions provided in these files decompress the input file in larger chunks
 * (using an internal buffer) and allow the user to easily get the lines one by
 * one, in order, by calling `gz_linebuffer_gets`. */

#ifndef GZ_READLINE_H
#define GZ_READLINE_H

#include <stddef.h>

/* Opens a file to read. At this time, only one file can be (globally) opened
 * for reading at a time. This function should not be called when another file
 * is open for reading.
 *
 * This function internally calls `gzopen`, so it works both on uncompressed
 * plaintext files and gzipped (.gz) files.
 *
 * Return value: 1 if file was opened successfully, 0 otherwise
 * */
int gz_linebuffer_open(const char *filename);

/* Close the file which is currently open.
 * */
void gz_linebuffer_close(void);

/* Read the next line of the currently open file and place it in `buffer`. This
 * function starts at the first line, and subsequent calls to it will yield the
 * second line, third line, etc. The newline at the end of the line is copied
 * into `buffer` as well (unless it's the last line and there is no trailing
 * newline), as well as a '\0' byte after it. If a line doesn't fit in the
 * buffer (or the newline and/or null byte don't fit), the function will trigger
 * a failed assert(). This will also happen if `gzread`, which is called
 * internally, returns a negative value (indicating an error happened).
 *
 * Return value: 1 if something was read and placed in `buffer`, 0 otherwise
 * (i.e. EOF was reached) */
int gz_linebuffer_gets(char *buffer, size_t bufsize);

#endif
