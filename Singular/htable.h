#ifndef HTABLE_H
#define HTABLE_H
/****************************************
*  Computer Algebra System SINGULAR     *
****************************************/

#include "kernel/mod2.h"

#include "Singular/subexpr.h"

struct stelem;
typedef struct stelem *   telem;

struct stelem
{
  telem next;
  char* key;
  sleftv val;
  uint32_t hash;
};

struct stablerec
{
  telem *t;
  int max;
  int ref;
};
#define TABLE_MAX 1000

stablerec* t_createTable(int s);
void t_destroyTable(stablerec* t);
stablerec* copyTable(stablerec* t);
char* stringTable(stablerec* t);
void htable_Print(stablerec *d);
/// find the entry to key s
telem t_findTable(stablerec* t,const char *s);
/// find the data to key s
leftv t_findTabelVal(stablerec* t,const char *s);
/// add a new entry (key s, data v) to table t
void t_addTable(stablerec* t,char *s, leftv v);
#endif
