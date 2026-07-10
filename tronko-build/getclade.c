/*
 *	getclade.c
 */
#include "getclade.h"
#include <string.h>

int getcladeArr(FILE *tre, struct masterArr *m, int max_nodename) {
    int n1, n2, n3;
    char ch;
    n1 = -1;
    n2 = -1;
    n3 = -1;
    do {
        if (specsearchArr(tre, m, max_nodename) == 1) {
            return tip + (m->numspec - 1);
        }
        ch = fgetc(tre);
        if (ch == ',') {
            comma++;
        }
        if (ch == ')') {
            if ((ch = fgetc(tre)) != ':') {
                ungetc(ch, tre);
            } else {
                do {
                    ch = (fgetc(tre));
                } while ((ch == '\n') || (ch == ' '));
                ungetc(ch, tre);
                fscanf(tre, "%lf", &m->tree[n3 - 1].bl);
            }
            return n3;
        }
        if (ch == '(') {
            n3 = getnodenumbArr(tre);
            n1 = getcladeArr(tre, m, max_nodename);
            n2 = getcladeArr(tre, m, max_nodename);
            if (n1 != -1 && n2 != -1 && n3 != -1) {
                linknodesArr(n1 - 1, n2 - 1, n3 - 1, m);
            }
            if (n1 == -1 || n2 == -1 || n3 == -1) {
                ungetc(ch, tre);
            }
        }
    } while (ch != ';');
    printf("Something went wrong (%s:%d)! Aborting...\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
}

void linknodesArr(int i, int j, int node, struct masterArr *m) {
    m->tree[node].up[0] = j;
    m->tree[node].up[1] = i;
    m->tree[i].down = node;
    m->tree[j].down = node;
}

int getnodenumbArr(FILE *tre) {
    char c;
    int i, j = 0;
    fpos_t position;
    i = 0;
    fgetpos(tre, &position);
    do {
        c = fgetc(tre);
        if (c == ',') i++;
        if (c == '(') j = j - 1;
        if (c == ')') j++;
    } while ((j < 0) && (c != EOF));
    fsetpos(tre, &position);
    return (i + comma + 1);
}

int specsearchArr(FILE *tre, struct masterArr *m, int max_node_name) {
    char ch;
    int i = 0;
    char specname[max_node_name + 1];
    ch = fgetc(tre);
    if ((ch != ')') && (ch != '(') && (ch != ',') && (ch != ' ') &&
        (ch != '\t') && (ch != '\n') && (ch != EOF)) {
        ungetc(ch, tre);
        while ((ch = fgetc(tre)) != ':' && (i < max_node_name)) {
            if (isalpha(ch) || isdigit(ch) || ch == '.' || ch == '_' ||
                ch == '/' || ch == '-' || ch == '|' || ch == '?' || ch == '*' ||
                ch == '&' || ch == '+' || ch == '#' || ch == '\'' ||
                ch == '=' || ch == ',' || ch == ']' || ch == '[' || ch == '>' ||
                ch == '<') {
                specname[i] = ch;
                i++;
            }
        }
        specname[i] = '\0';
        int tmp1 = 0;
        for (tmp1 = 0; tmp1 < m->numspec; tmp1++) {
            if (!strcmp(m->names[tmp1], specname)) {
                tip = tmp1 + 1;
            }
        }
        strcpy(m->tree[tip + m->numspec - 2].name, specname);
        m->tree[tip + m->numspec - 2].up[0] = -1;
        m->tree[tip + m->numspec - 2].up[1] = -1;
        fscanf(tre, "%lf", &m->tree[tip + m->numspec - 2].bl);
        return 1;
    } else {
        ungetc(ch, tre);
        return 0;
    }
}
