/** \brief
 * C language version of MATLAB Wavelet Toolbox
 * the essential wave is sym8
 * \return
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "List.h"
#include "Vector2.h"
#define DATA_LEN 256
typedef unsigned char           uint8_t;
typedef unsigned short          uint16_t;

enum ConvolutionType
{
    kFullConvolution = 0,
    kSameConvolution,
    kValidConvolution
};
//function declaration
void Wavedec(List_Exten * x, uint16_t num_levels);
void Dwt(List_Exten * x,List_Exten * a,List_Exten * d);
List_Exten* List_Conv(List_Exten *v,const List_Exten* filter,enum ConvolutionType convolutiontype);
List_Exten* Wrcoef(List_Exten* l,Vector2*c,uint16_t level);
List_Exten* Convdown(List_Exten* x,List_Exten* f);
List_Exten* Upconv(List_Exten* x,List_Exten* f, uint16_t s);
List_Exten* Dyadup(List_Exten* x);
List_Exten* Wextend(List_Exten* x, uint16_t lenEXT);
List_Exten* Dwta(List_Exten * x);
List_Exten* Dwtd(List_Exten * x);
//--------------Global variable declaration----------
//sym8's filter parameter
const float Lo_D[16] =
{
    -0.0033824159510061256,
        -0.00054213233179114812,
        0.031695087811492981,
        0.0076074873249176054,
        -0.14329423835080971,
        -0.061273359067658524,
        0.48135965125837221,
        0.77718575170052351,
        0.3644418948353314,
        -0.051945838107709037,
        -0.027219029917056003,
        0.049137179673607506,
        0.0038087520138906151,
        -0.014952258337048231,
        -0.0003029205147213668,
        0.0018899503327594609
    };
const float Hi_D[16] =
{
    -0.0018899503327594609,
        -0.0003029205147213668,
        0.014952258337048231,
        0.0038087520138906151,
        -0.049137179673607506,
        -0.027219029917056003,
        0.051945838107709037,
        0.3644418948353314,
        -0.77718575170052351,
        0.48135965125837221,
        0.061273359067658524,
        -0.14329423835080971,
                -0.0076074873249176054,
        0.031695087811492981,
        0.00054213233179114812,
        -0.0033824159510061256
    };
const float Lo_R[16] =
{
    0.0018899503327594609,
    -0.0003029205147213668,
    -0.014952258337048231,
    0.0038087520138906151,
    0.049137179673607506,
    -0.027219029917056003,
    -0.051945838107709037,
    0.3644418948353314,
    0.77718575170052351,
    0.48135965125837221,
    -0.061273359067658524,
    -0.14329423835080971,
    0.0076074873249176054,
    0.031695087811492981,
    -0.00054213233179114812,
    -0.0033824159510061256
};
const float Hi_R[16] =
{
    -0.0033824159510061256,
        0.00054213233179114812,
        0.031695087811492981,
        -0.0076074873249176054,
        -0.14329423835080971,
        0.061273359067658524,
        0.48135965125837221,
        -0.77718575170052351,
        0.3644418948353314,
        0.051945838107709037,
        -0.027219029917056003,
        -0.049137179673607506,
        0.0038087520138906151,
        0.014952258337048231,
        -0.0003029205147213668,
        -0.0018899503327594609
    };
//Global variable declaration
const uint16_t xdata = 8;
Vector2 *C;
List_Exten* L;
List_Exten * Hi_R_;
List_Exten * Lo_R_;
List_Exten * Hi_D_;
List_Exten * Lo_D_;
int main()
{
    //Converts the filter parameter array to a custom List structure
    Lo_D_ = (List_Exten*)malloc(sizeof(List_Exten));
    List_Init(Lo_D_,16);
    Hi_D_ = (List_Exten*)malloc(sizeof(List_Exten));
    List_Init(Hi_D_,16);
    Lo_R_ = (List_Exten*)malloc(sizeof(List_Exten));
    List_Init(Lo_R_,16);
    Hi_R_ = ( List_Exten*)malloc(sizeof(List_Exten));
    List_Init(Hi_R_,16);
    for(uint16_t i =0; i<16; ++i)
    {
        List_Append(Lo_D_,Lo_D[i]);
        List_Append(Hi_D_,Hi_D[i]);
        List_Append(Lo_R_,Lo_R[i]);
        List_Append(Hi_R_,Hi_R[i]);
    }
    //Declaration of c and l according to the Matlab version
    L = (List_Exten*)malloc(sizeof(List_Exten));
    List_Init(L,xdata+1);
    C = (Vector2*)malloc(sizeof( Vector2));
    Vector2_Init(C,xdata+1);

    //Example data , for sym8, the data in principle must be greater than 256 elements to be meaningful.
    float data[DATA_LEN] = {10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 21, 22};
    List_Exten *listpressure = (List_Exten*)malloc(sizeof(List_Exten));
    List_Init(listpressure,DATA_LEN);

    for(uint16_t i =0; i<DATA_LEN; ++i)
    {
        List_Append(listpressure,data[i]);
    }

    //The results are automatically assigned to the global variables C and L
    Wavedec(listpressure,xdata);
    free(listpressure);
    //After the wavelet decomposition C is directly C in MATLAB, and L is inserted after the original data length  in MATLAB
    List_Append(L,DATA_LEN);

    //Wavelet reconstruction, the selected reconstruction level is the seventh layer
    List_Exten* result = Wrcoef(L,C,7);
    List_Show(result);

    free(C);
    free(L);
    free(Lo_D_);
    free(Hi_D_);
    free(Lo_R_);
    free(Hi_R_);
    return 0;
}
/** \brief
 *  Wavelet Decompositoin
 * \param x is the original input data
 * \param num_levels is the level to decompse
 * \return
 *
 */
void Wavedec(List_Exten * x, uint16_t num_levels)
{
    if(x->n < (1 << num_levels))
    {
        printf("std::invalid_argument(Wavelet::Wavedec)");
    }
    //a
    List_Exten* a = (List_Exten*)malloc(sizeof(List_Exten));
    List_Init(a,x->n);
    for (uint16_t i= 0; i < x->n; ++i)
    {
        List_Append(a,x->A[i]);

    }
    Vector2 *data = (Vector2*)malloc(sizeof(Vector2));
    Vector2_Init(data,num_levels+1);
    List_Exten *len = (List_Exten*)malloc(sizeof(List_Exten));
    List_Init(len,num_levels+1);

    for(uint16_t k = 0 ; k < num_levels ; ++k )
    {

        List_Exten * ak = (List_Exten*)malloc(sizeof(List_Exten));
        List_Exten * dk = (List_Exten*)malloc(sizeof(List_Exten));

        ak=Dwta(a);
        dk=Dwtd(a);

        free(a->A);
        free(a);
        a = (List_Exten*)malloc(sizeof(List_Exten));
        List_Init(a,ak->n);
        for (uint16_t i = 0; i < ak->n; ++i)
        {
            List_Append(a,ak->A[i]);
        }
        Vector2_Append(data,dk);
        List_Append(len,dk->n);
    }
    Vector2_Append(data,a);
    free(a);
    List_Append(len,a->n);

    for(uint16_t k=data->n; k>=1; --k)
    {
        Vector2_Append(C,Vector2_GetList(data,k-1));
    }
    for(uint16_t k=len->n; k>=1; --k)
    {
        List_Append(L,len->A[k-1]);
    }
    free(data);
    free(len);
}
List_Exten* Dwta(List_Exten * x)
{
    List_Exten* data = (List_Exten*)malloc(sizeof(List_Exten));
    List_Init(data,4);
    data = Convdown(x,Lo_D_);
    return data;
}

List_Exten* Dwtd(List_Exten * x)
{
    List_Exten* data = (List_Exten*)malloc(sizeof(List_Exten));
    List_Init(data,4);
    data = Convdown(x,Hi_D_);
    return data;
}

List_Exten* Convdown(List_Exten* x, List_Exten* f)
{
    uint16_t last = (uint16_t)(ceil((float)(x->n+1) / 2.0) * 2.0);
    List_Exten* z = Wextend(x, xdata);
    z = List_Conv(z, f, kValidConvolution);
    List_Exten* data=(List_Exten*)malloc(sizeof(List_Exten));
    List_Init(data,last/2);
    for (uint16_t k = 1; k <= last; k += 2)
    {
        List_Append(data,z->A[k]);
    }
    free(z->A);
    free(z);
    return data;
}
/** \brief
 * Wavelet reconstruction function
 * \param l is the result of Wavelet Decomposition
 * \param c is the result of Wavelet Decomposition, the length of each element is storaged in l.
 * \param level is the level in order to reconstruct
 * \return
 *
 */
List_Exten* Wrcoef(List_Exten* l,Vector2*c,uint16_t level)
{
    uint16_t iMin = l->n-level;
    List_Exten* data = Vector2_GetList(c,iMin-1);
    data= Upconv(data,Hi_R_,l->A[iMin]);
    for(uint16_t k=1; k<level ; ++k)
    {
        data = Upconv(data,Lo_R_,l->A[iMin+k]);
    }
    return data;
}

List_Exten* List_Conv(List_Exten *v,const List_Exten* filter,enum ConvolutionType convolutiontype)
{
    uint16_t vsize = v->n;
    uint16_t fsize = filter->n;

    if (fsize > vsize)
    {
        printf("throw std::invalid_argument(conv)\n");
    }
    uint16_t mini = 0;
    uint16_t maxi = vsize + fsize - 1;
    if (convolutiontype == kSameConvolution)
    {
        mini = fsize / 2;
        maxi = mini + vsize;
    }
    else if (convolutiontype == kValidConvolution)
    {
        mini = fsize - 1;
        maxi = vsize;
    }
    List_Exten* data=(List_Exten*)malloc(sizeof(List_Exten));
    List_Init(data,maxi - mini);
    for (uint16_t i = mini; i < maxi; i++)
    {
        uint16_t minj = (i < fsize - 1) ? 0 : i - fsize + 1;
        uint16_t maxj = (i >= vsize - 1) ? vsize : i + 1;
        float value =0;
        for (uint16_t j = minj; j < maxj; j++)
        {
            value += v->A[j] * filter->A[i - j];
        }
        List_Append(data,value);
    }

    return data;
}


List_Exten* Upconv(List_Exten* x,List_Exten* f,uint16_t s)
{
    uint16_t lf = f->n;
    List_Exten* data = Dyadup(x);
    data = Wextend(data, lf / 2);
    data = List_Conv(data, f, kFullConvolution);
    List_erase(data,0,lf - 1);
    List_erase(data,s,data->n);
    return data;
}

List_Exten* Dyadup(List_Exten* x)
{
    List_Exten *data=(List_Exten*)malloc(sizeof(List_Exten));
    List_Init(data,2 * x->n);
    for (uint16_t k = 0; k < x->n; k++)
    {
        List_Append(data,x->A[k]);
        List_Append(data,0);
    }
    return data;
}

List_Exten* Wextend(List_Exten* x, uint16_t lenEXT)
{
    List_Exten* data = (List_Exten*)malloc(sizeof(List_Exten));
    List_Init(data,x->n);
    for (uint16_t i = 0; i < x->n; ++i)
    {
        List_Append(data,x->A[i]);

    }
    if ((x->n )% 2 == 1)
    {
        List_Append(data,data->A[data->n-1]);//temp.push_back(temp.back());
    }
    uint16_t rep = lenEXT / (data->n);
    uint16_t len = lenEXT % (data->n);
    List_Exten* y = (List_Exten*)malloc(sizeof(List_Exten));
    List_Init(y,2 * (len + rep * (data->n)) + data->n);
    for(uint16_t i=data->n-len; i<data->n; i++)
    {
        List_Append(y,data->A[i]);
    }

    for(uint16_t k = 0; k < 2 * rep + 1; k++)
    {
        for(uint16_t i=0; i<data->n; i++)
        {
            List_Append(y,data->A[i]);
        }
    }
    for(uint16_t i=0; i<len; i++)
        for(uint16_t i=0; i<len; i++)
        {
            List_Append(y,data->A[i]);
        }
    free(data->A);
    free(data);
    return y;
}
