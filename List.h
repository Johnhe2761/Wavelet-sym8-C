#ifndef _h_List_
#define _h_List_

typedef unsigned char           uint8_t;
typedef unsigned short          uint16_t;
typedef struct List_Exten
{
    float *A;
    uint16_t n;
    uint16_t c; // capacity
} List_Exten;
void List_Init(List_Exten* List,uint16_t capacity);
void List_Append(List_Exten* List, float value);
void List_Show(List_Exten* List);
List_Exten* List_Add(List_Exten* L1,List_Exten* L2);
void List_erase(List_Exten *v,uint16_t vbegin,uint16_t vend);
#endif
