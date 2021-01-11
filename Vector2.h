#ifndef _h_Vector2_
#define _h_Vector2_
#include "List.h"

typedef unsigned char           uint8_t;
typedef unsigned short          uint16_t;
typedef struct Vector2
{
    struct List_Exten * vec;
    uint16_t n;
    uint16_t c;
} Vector2;
void Vector2_Init(Vector2* vec2,uint16_t capacity);
void Vector2_Append(Vector2* vec2,List_Exten* value);
void Vector2_Show(Vector2* vec2);
List_Exten* Vector2_GetList(Vector2 * v, uint16_t x);

#endif

