
#include "List.h"
/** \brief
 *  initialization of the struct
 * \param data is a initialized struct
 * \param capacity is the initial size of the struct
 * \return
 *
 */
void List_Init(List_Exten* data,uint16_t capacity)
{
    data->c = capacity;
    data->n = 0;
    data->A = (float*)malloc(sizeof(float)*data->c);
}
/** \brief
 *  Add element into the struct
 * \param data is the target to add
 * \param value is the element, if the capacity is not enough, it
 *         will expand based on the initial size automatically
 * \return
 *
 */
void List_Append(List_Exten* data, float value)
{
    if (data->c == data->n)
    {
        data->c = 2 * data->c;
        float* B = (float*)malloc(sizeof(float)*data->c);
        for (uint16_t i = 0; i < data->n; ++i)
        {
            B[i] = data->A[i];
        }
        //printf("extented!\n");
        free(data->A);
        data->A = B;
    }
    data->A[data->n++] = value;
}

void List_Show(List_Exten* data)
{
    for (uint16_t i = 0; i < data->n; ++i)
    {
        printf("%lf\n",data->A[i]);
    }
}
/** \brief
 *  merge two List structs into one;
 * \return
 *
 */

List_Exten* List_Add(List_Exten* L1,List_Exten* L2)
{
    if(L1->n!=L2->n)
    {
        printf("not the same size!\n");
        return 0;
    }
    List_Exten* data = (List_Exten*)malloc(sizeof(List_Exten));
    List_Init(data,L1->n);
    for(uint16_t i=0; i<L1->n; i++)
    {
        List_Append(data,L1->A[i]+L2->A[i]);
    }
    return data;
}
/** \brief
 *  remove a piece of data from vbegin to vend
 * \return
 *
 */

void List_erase(List_Exten *v,uint16_t vbegin,uint16_t vend)
{
    v->n=v->n-vend+vbegin;
    float* data = (float*)malloc(sizeof(float)*v->n);
    if(vbegin==0)
    {
        for (uint16_t i = 0; i < v->n; ++i)
        {
            data[i] = v->A[vend+i];
        }
    }
    else
    {
        for (uint16_t i = 0; i < v->n; ++i)
        {
            data[i] = v->A[i];
        }
    }
    free(v->A);
    v->A = data;
}
