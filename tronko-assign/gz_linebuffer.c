#include "gz_linebuffer.h"
#include <zlib.h>
#include <string.h>
#include <assert.h>

#define GZ_INTERNAL_BUF_SIZE (1 << 20) // = 1 MiB

static gzFile gz_file = Z_NULL;
static unsigned char gz_inbuf[GZ_INTERNAL_BUF_SIZE];
static size_t gz_inbuf_pos = 0, gz_inbuf_len = 0;

int gz_linebuffer_open(const char *filename) {
    assert(Z_NULL == gz_file);
    gz_file = gzopen(filename, "r");
    gz_inbuf_pos = gz_inbuf_len = 0;
    return Z_NULL != gz_file;
}

void gz_linebuffer_close(void) {
    assert(gz_file != Z_NULL);
    gzclose(gz_file);
    gz_file = Z_NULL;
}

int gz_linebuffer_gets(char *buffer, size_t bufsize) {
    size_t out_pos = 0;

    while (1) {
        if (gz_inbuf_pos == gz_inbuf_len) {
            int n = gzread(gz_file, gz_inbuf, GZ_INTERNAL_BUF_SIZE);
            assert(n >= 0 && "gz_readline: gzread error");
            if (n == 0) {
                if (out_pos) buffer[out_pos] = '\0';
                return out_pos > 0;
            }
            gz_inbuf_len = (size_t)n;
            gz_inbuf_pos = 0;
        }

        unsigned char *start = gz_inbuf + gz_inbuf_pos;
        unsigned char *nl = memchr(start, '\n', gz_inbuf_len - gz_inbuf_pos);
        size_t chunk =
            nl ? (size_t)(nl - start) + 1 : gz_inbuf_len - gz_inbuf_pos;

        assert(out_pos + chunk < bufsize &&
               "gz_readline: line too long for buffer");
        memcpy(buffer + out_pos, start, chunk);
        out_pos += chunk;
        gz_inbuf_pos += chunk;

        if (nl) {
            buffer[out_pos] = '\0';
            return 1;
        }
    }
}
