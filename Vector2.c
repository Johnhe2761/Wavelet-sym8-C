
#include "Vector2.h"
/** \brief
 *  initialization of the struct
 * \param data is a initialized struct
 * \param capacity is the initial size of the struct
 * \return
 *
 */

void Vector2_Init(Vector2* data,uint16_t capacity)
{
    data->c = capacity;
    data->n = 0;
    data->vec = (List_Exten *)malloc(sizeof(List_Exten)*(data->c));
}
/** \brief
 *  Add element into the struct
 * \param data is the target to add
 * \param value is the element, if the capacity is not enough, it
 *         will expand based on the initial size automatically
 * \return
 *
 */
void Vector2_Append(Vector2* data,List_Exten* value)
{
    if (data->c == data->n)
    {
        data->c = 2 * data->c;
        List_Exten* B = (List_Exten*)malloc(sizeof(List_Exten)*data->c);
        for (uint16_t i = 0; i < data->n; ++i)
        {
            B[i] = data->vec[i];
        }
        //printf("Vecotor2_Append()\n");
        free(data->vec);
        data->vec = B;
    }
    data->vec[data->n++] = *value;
}
void Vector2_Show(Vector2* data)
{
    for (uint16_t i = 0; i < data->n; ++i)
    {
        List_Show(&data->vec[i]);
    }
}
/** \brief
 *  Get the element according to the index given
 * \param v is the target struct
 * \param x is the index
 * \return
 *
 */

List_Exten* Vector2_GetList(Vector2 * v, uint16_t x)
{
    List_Exten * data = (List_Exten*)malloc(sizeof(List_Exten));
    List_Init(data,4);
    for(uint16_t k = 0 ; k<v->vec[x].n; ++k)
    {
        List_Append(data,v->vec[x].A[k]);
    }
    return data;
}

