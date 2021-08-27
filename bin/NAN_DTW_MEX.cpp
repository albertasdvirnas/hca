/***********************************************************************/
/************************* DISCLAIMER **********************************/
/***********************************************************************/
/** This UCR Suite software is copyright protected � 2012 by          **/
/** Thanawin Rakthanmanon, Bilson Campana, Abdullah Mueen,            **/
/** Gustavo Batista and Eamonn Keogh.                                 **/
/**                                                                   **/
/** Unless stated otherwise, all software is provided free of charge. **/
/** As well, all software is provided on an "as is" basis without     **/
/** warranty of any kind, express or implied. Under no circumstances  **/
/** and under no legal theory, whether in tort, contract,or otherwise,**/
/** shall Thanawin Rakthanmanon, Bilson Campana, Abdullah Mueen,      **/
/** Gustavo Batista, or Eamonn Keogh be liable to you or to any other **/
/** person for any indirect, special, incidental, or consequential    **/
/** damages of any character including, without limitation, damages   **/
/** for loss of goodwill, work stoppage, computer failure or          **/
/** malfunction, or for any and all other damages or losses.          **/
/**                                                                   **/
/** If you do not agree with these terms, then you you are advised to **/
/** not use this software.                                            **/
/***********************************************************************/
/***********************************************************************/

// this version includes support for nan's

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <iostream>
#include "mex.h"


#define min(x,y) ((x)<(y)?(x):(y))
#define max(x,y) ((x)>(y)?(x):(y))
#define dist(x,y) ((x-y)*(x-y))

#define INF 1e20       //Pseudo Infitinte number for this code

using namespace std;

/// Data structure for sorting the query
typedef struct Index
    {   double value;
        int    index;
    } Index;

/// Data structure (circular array) for finding minimum and maximum for LB_Keogh envolop
struct deque
{   int *dq;
    int size,capacity;
    int f,r;
};


/// Sorting function for the query, sort by abs(z_norm(q[i])) from high to low
int comp(const void *a, const void* b)
{   Index* x = (Index*)a;
    Index* y = (Index*)b;
    return abs(y->value) - abs(x->value);   // high to low
}

/// Initial the queue at the begining step of envelop calculation
void init(deque *d, int capacity)
{
    d->capacity = capacity;
    d->size = 0;
    d->dq = (int *) malloc(sizeof(int)*d->capacity);
    d->f = 0;
    d->r = d->capacity-1;
}

/// Destroy the queue
void destroy(deque *d)
{
    free(d->dq);
}

/// Insert to the queue at the back
void push_back(struct deque *d, int v)
{
    d->dq[d->r] = v;
    d->r--;
    if (d->r < 0)
        d->r = d->capacity-1;
    d->size++;
}

/// Delete the current (front) element from queue
void pop_front(struct deque *d)
{
    d->f--;
    if (d->f < 0)
        d->f = d->capacity-1;
    d->size--;
}

/// Delete the last element from queue
void pop_back(struct deque *d)
{
    d->r = (d->r+1)%d->capacity;
    d->size--;
}

/// Get the value at the current position of the circular queue
int front(struct deque *d)
{
    int aux = d->f - 1;

    if (aux < 0)
        aux = d->capacity-1;
    return d->dq[aux];
}

/// Get the value at the last position of the circular queueint back(struct deque *d)
int back(struct deque *d)
{
    int aux = (d->r+1)%d->capacity;
    return d->dq[aux];
}

/// Check whether or not the queue is empty
int empty(struct deque *d)
{
    return d->size == 0;
}

/// Finding the envelop of min and max value for LB_Keogh
/// Implementation idea is intoruduced by Danial Lemire in his paper
/// "Faster Retrieval with a Two-Pass Dynamic-Time-Warping Lower Bound", Pattern Recognition 42(9), 2009.
void lower_upper_lemire(double *t, int len, int r, double *l, double *u)
{
    struct deque du, dl;

    init(&du, 2*r+2);
    init(&dl, 2*r+2);

    push_back(&du, 0);
    push_back(&dl, 0);

    for (int i = 1; i < len; i++)
    {
        if (i > r)
        {
            u[i-r-1] = t[front(&du)];
            l[i-r-1] = t[front(&dl)];
        }
        if (t[i] > t[i-1])
        {
            pop_back(&du);
            while (!empty(&du) && t[i] > t[back(&du)])
                pop_back(&du);
        }
        else
        {
            pop_back(&dl);
            while (!empty(&dl) && t[i] < t[back(&dl)])
                pop_back(&dl);
        }
        push_back(&du, i);
        push_back(&dl, i);
        if (i == 2 * r + 1 + front(&du))
            pop_front(&du);
        else if (i == 2 * r + 1 + front(&dl))
            pop_front(&dl);
    }
    for (int i = len; i < len+r+1; i++)
    {
        u[i-r-1] = t[front(&du)];
        l[i-r-1] = t[front(&dl)];
        if (i-front(&du) >= 2 * r + 1)
            pop_front(&du);
        if (i-front(&dl) >= 2 * r + 1)
            pop_front(&dl);
    }
    destroy(&du);
    destroy(&dl);
}

/// Calculate quick lower bound
/// Usually, LB_Kim take time O(m) for finding top,bottom,fist and last.
/// However, because of z-normalization the top and bottom cannot give siginifant benefits.
/// And using the first and last points can be computed in constant time.
/// The prunning power of LB_Kim is non-trivial, especially when the query is not long, say in length 128.
double lb_kim_hierarchy(double *t, double *q, int j, int len, double mean, double std, double bsf = INF)
{
    /// 1 point at front and back
    double d, lb;
    double x0 = (t[j] - mean) / std;
    double y0 = (t[(len-1+j)] - mean) / std;
    lb = dist(x0,q[0]) + dist(y0,q[len-1]);
    if (lb >= bsf)   return lb;

    /// 2 points at front
    double x1 = (t[(j+1)] - mean) / std;
    d = min(dist(x1,q[0]), dist(x0,q[1]));
    d = min(d, dist(x1,q[1]));
    lb += d;
    if (lb >= bsf)   return lb;

    /// 2 points at back
    double y1 = (t[(len-2+j)] - mean) / std;
    d = min(dist(y1,q[len-1]), dist(y0, q[len-2]) );
    d = min(d, dist(y1,q[len-2]));
    lb += d;
    if (lb >= bsf)   return lb;

    /// 3 points at front
    double x2 = (t[(j+2)] - mean) / std;
    d = min(dist(x0,q[2]), dist(x1, q[2]));
    d = min(d, dist(x2,q[2]));
    d = min(d, dist(x2,q[1]));
    d = min(d, dist(x2,q[0]));
    lb += d;
    if (lb >= bsf)   return lb;

    /// 3 points at back
    double y2 = (t[(len-3+j)] - mean) / std;
    d = min(dist(y0,q[len-3]), dist(y1, q[len-3]));
    d = min(d, dist(y2,q[len-3]));
    d = min(d, dist(y2,q[len-2]));
    d = min(d, dist(y2,q[len-1]));
    lb += d;

    return lb;
}

/// LB_Keogh 1: Create Envelop for the query
/// Note that because the query is known, envelop can be created once at the begenining.
///
/// Variable Explanation,
/// order : sorted indices for the query.
/// uo, lo: upper and lower envelops for the query, which already sorted.
/// t     : a circular array keeping the current data.
/// j     : index of the starting location in t
/// cb    : (output) current bound at each position. It will be used later for early abandoning in DTW.
double lb_keogh_cumulative(int* order, double *t, double *uo, double *lo, double *cb, int j, int len, double mean, double std, double best_so_far = INF)
{
    double lb = 0;
    double x, d;

    for (int i = 0; i < len && lb < best_so_far; i++)
    {
        x = (t[(order[i]+j)] - mean) / std;
        d = 0;
        if (x > uo[i])
            d = dist(x,uo[i]);
        else if(x < lo[i])
            d = dist(x,lo[i]);
        lb += d;
        cb[order[i]] = d;
    }
    return lb;
}

/// LB_Keogh 2: Create Envelop for the data
/// Note that the envelops have been created (in main function) when each data point has been read.
///
/// Variable Explanation,
/// tz: Z-normalized data
/// qo: sorted query
/// cb: (output) current bound at each position. Used later for early abandoning in DTW.
/// l,u: lower and upper envelop of the current data
double lb_keogh_data_cumulative(int* order, double *tz, double *qo, double *cb, double *l, double *u, int len, double mean, double std, double best_so_far = INF)
{
    double lb = 0;
    double uu,ll,d;

    for (int i = 0; i < len && lb < best_so_far; i++)
    {
        uu = (u[order[i]]-mean)/std;
        ll = (l[order[i]]-mean)/std;
        d = 0;
        if (qo[i] > uu)
            d = dist(qo[i], uu);
        else
        {   if(qo[i] < ll)
            d = dist(qo[i], ll);
        }
        lb += d;
        cb[order[i]] = d;
    }
    return lb;
}

/// Calculate Dynamic Time Wrapping distance
/// A,B: data and query, respectively
/// cb : cummulative bound used for early abandoning
/// r  : size of Sakoe-Chiba warpping band
double dtw(double* A, double* B, double *cb, int m, int r, double bsf = INF)
{

    double *cost;
    double *cost_prev;
    double *cost_tmp;
    int i,j,k;
    double x,y,z,min_cost;

    /// Instead of using matrix of size O(m^2) or O(mr), we will reuse two array of size O(r).
    cost = (double*)malloc(sizeof(double)*(2*r+1));
    for(k=0; k<2*r+1; k++)    cost[k]=INF;

    cost_prev = (double*)malloc(sizeof(double)*(2*r+1));
    for(k=0; k<2*r+1; k++)    cost_prev[k]=INF;

    for (i=0; i<m; i++)
    {
        k = max(0,r-i);
        min_cost = INF;

        for(j=max(0,i-r); j<=min(m-1,i+r); j++, k++)
        {
            /// Initialize all row and column
            if ((i==0)&&(j==0))
            {
                cost[k]=dist(A[0],B[0]);
                min_cost = cost[k];
                continue;
            }

            if ((j-1<0)||(k-1<0))     y = INF;
            else                      y = cost[k-1];
            if ((i-1<0)||(k+1>2*r))   x = INF;
            else                      x = cost_prev[k+1];
            if ((i-1<0)||(j-1<0))     z = INF;
            else                      z = cost_prev[k];

            /// Classic DTW calculation
            cost[k] = min( min( x, y) , z) + dist(A[i],B[j]);

            /// Find minimum cost in row for early abandoning (possibly to use column instead of row).
            if (cost[k] < min_cost)
            {   min_cost = cost[k];
            }
        }

        /// We can abandon early if the current cummulative distace with lower bound together are larger than bsf
        if (i+r < m-1 && min_cost + cb[i+r+1] >= bsf)
        {   free(cost);
            free(cost_prev);
            return min_cost + cb[i+r+1];
        }

        /// Move current array to previous array.
        cost_tmp = cost;
        cost = cost_prev;
        cost_prev = cost_tmp;
    }
    k--;

    /// the DTW distance is in the last cell in the matrix of size O(m^2) or at the middle of our array.
    double final_dtw = cost_prev[k];
    free(cost);
    free(cost_prev);
    return final_dtw;
}
//
// double TwoSum(double a,double b, double *l, double *u)
// {
//     double x = a + b;
//     double z = x - a;
//     double y = (a - (x - z))+(b-z);
//     return x,y;
// }
/// Print function for debugging
void printArray(double *x, int len)
{   for(int i=0; i<len; i++)
        printf(" %6.2lf",x[i]);
    printf("\n");
}

/// If expected error happens, teminated the program.
void error(int id)
{
    if(id==1)
        printf("ERROR : Memory can't be allocated!!!\n\n");
    else if ( id == 2 )
        printf("ERROR : File not Found!!!\n\n");
    else if ( id == 3 )
        printf("ERROR : Can't create Output File!!!\n\n");
    else if ( id == 4 )
    {
        printf("ERROR : Invalid Number of Arguments!!!\n");
        ///printf("Command Usage:  UCR_DTW.exe  data-file  query-file   m   R\n\n");
        ///printf("For example  :  UCR_DTW.exe  data.txt   query.txt   128  0.05\n");
    }
    ///exit(1);
}

class mystream : public std::streambuf
{
protected:
virtual std::streamsize xsputn(const char *s, std::streamsize n) { mexPrintf("%.*s", n, s); return n; }
virtual int overflow(int c=EOF) { if (c != EOF) { mexPrintf("%.1s", &c); } return 1; }
};
class scoped_redirect_cout
{
public:
  scoped_redirect_cout() { old_buf = std::cout.rdbuf(); std::cout.rdbuf(&mout); }
  ~scoped_redirect_cout() { std::cout.rdbuf(old_buf); }
private:
  mystream mout;
  std::streambuf *old_buf;
};




/// Main Function
/// int main(  int argc , char *argv[] )
void cpp_dtw(double *data, double *q, double *b, int n,int m, double R, double *y, double *s)
{

    scoped_redirect_cout mycout_redirect;

    double bsf;      /// best-so-far
    double *t;       /// data array and query array
    int *order;      ///new order of the query
    double *u, *l, *qo, *uo, *lo,*tz,*cb, *cb1, *cb2,*u_d, *l_d;

    // outputs..
    // y = 0;

    *s = 0;


    /// name the parameters
    double d;
    long long i , j;
    double ex , ex2 , mean, std;
    int r = -1;
    long long loc = 0;
    double t1,t2;
    int kim = 0,keogh = 0, keogh2 = 0; int bm=0;
    double dist=0, lb_kim=0, lb_k=0, lb_k2=0;
    double *buffer,*u_buff, *l_buff;
    double *bufferb,*tb;
    Index *Q_tmp;

    /// For every EPOCH points, all cummulative values, such as ex (sum), ex2 (sum square), will be restarted for reducing the floating point error.
    int EPOCH = 100000;

    /// Sakoe-Chiba band
    if (R<=1)
            r = floor(R*m);
	else
            r = floor(R);
    // printf("%d\n",r);



    /// start the clock
    t1 = clock();

// initializations
{
    qo = (double *)malloc(sizeof(double)*m);
    uo = (double *)malloc(sizeof(double)*m);
    lo = (double *)malloc(sizeof(double)*m);
    order = (int *)malloc(sizeof(int)*m);
    Q_tmp = (Index *)malloc(sizeof(Index)*m);
    u = (double *)malloc(sizeof(double)*m);
    l = (double *)malloc(sizeof(double)*m);
    cb = (double *)malloc(sizeof(double)*m);
    cb1 = (double *)malloc(sizeof(double)*m);
    cb2 = (double *)malloc(sizeof(double)*m);
    u_d = (double *)malloc(sizeof(double)*m);
    l_d = (double *)malloc(sizeof(double)*m);
    t = (double *)malloc(sizeof(double)*m*2);
    tb = (double *)malloc(sizeof(double)*m*2);

    tz = (double *)malloc(sizeof(double)*m);
    buffer = (double *)malloc(sizeof(double)*EPOCH);
    bufferb = (double *)malloc(sizeof(double)*EPOCH);

    u_buff = (double *)malloc(sizeof(double)*EPOCH);
    l_buff = (double *)malloc(sizeof(double)*EPOCH);

    if( qo == NULL || uo == NULL || lo == NULL|| order == NULL|| Q_tmp == NULL || u == NULL || l == NULL || cb == NULL|| cb1 == NULL || cb2 == NULL ||u_d == NULL|| l_d == NULL || t == NULL ||  tz == NULL )
    {
        error(1);
        return;
    }
}




    /// Read query file
    bsf = INF; /// So this finds global solution
    j = 0;
    ex = ex2 = 0;

    // compute statistics for query /// alternatively use Error-free transformation here
    for( i = 0 ; i < m ; i++ )
    {
        ex += q[i];    /// if the numbers are very big, would have some numerical issues?
        ex2 += q[i]*q[i];
    }
    //Sum(a,b)

    /// Do z-normalize the query, keep in same array, q
    mean = ex/m;
    std = ex2/m;
    std = sqrt(std-mean*mean);
    for( i = 0 ; i < m ; i++ )
         q[i] = (q[i] - mean)/std;

    lower_upper_lemire(q, m, r, l, u);

    // printf("%f\n",l[0]);
    // printf("%f\n",u[0]);
    FILE *File2;
    File2 = fopen("lemire_l.txt", "w+");
    for(int k = 0 ; k < m ; k++ )
        fprintf(File2,"%f ",l[k]);
    fclose(File2);

    File2 = fopen("lemire_u.txt", "w+");
    for(int k = 0 ; k < m ; k++ )
        fprintf(File2,"%f ",u[k]);
    fclose(File2);



    // test purposes: save lower and upper lemire to txt files

    /// Sort the query one time by abs(z-norm(q[i]))
    for( i = 0; i<m; i++)
    {
        Q_tmp[i].value = q[i];
        Q_tmp[i].index = i;
    }
    qsort(Q_tmp, m, sizeof(Index),comp);

    /// also create another arrays for keeping sorted envelop
    for( i=0; i<m; i++)
    {   int o = Q_tmp[i].index;
        order[i] = o;
        qo[i] = q[o];
        uo[i] = u[o];
        lo[i] = l[o];
    }
    free(Q_tmp);


    /// Initial the cummulative lower bound / also consider uper bound maybe?
    for( i=0; i<m; i++)
    {   cb[i]=0;
        cb1[i]=0;
        cb2[i]=0;
    }

    i = 0;          /// current index of the data in current chunk of size EPOCH
    j = 0;          /// the starting index of the data in the circular array, t
    ex = ex2 = 0;
    bool done = false;
    int it=0, ep=0, k=0;
    long long I;    /// the starting index of the data in current chunk of size EPOCH
    double curbit = 0;
    /// Read first m-1 points
    // ep=0;
    // if (it==0)
    // {   for(k=0; k<m-1; k++)
    //         buffer[k] = data[k];
    //         bufferb[k] = b[k];
    //         // printf(" %4.3f \n",  bufferb[k]);
    //
    // }
    // printf("%f\n",buffer[0]);
    // printf(" %3d \n", b[99]);
    
//     FILE *File3;
//     File3 = fopen("out.txt", "w+");
    while(!done)
    {
        // printf("Iteration %d\n",it);
        /// Read first m-1 points
         // printf(" %3d \n",m);
        ep=0;
        if (it==0)
        {   for(k=0; k<m-1; k++){
                buffer[k] = data[k];
                bufferb[k] = b[k];
                // printf(" %3d \n", b[k]);
            }


        }
        else
        {   for(k=0; k<m-1; k++){
                buffer[k] = buffer[EPOCH-m+1+k];
                bufferb[k] = bufferb[EPOCH-m+1+k];
            }
        }

        /// Read buffer of size EPOCH or when all data has been read.
        ep=m-1; // create buffer.. since all data is already loaded, might not be necessary?
        while(ep<EPOCH && EPOCH*it+ep <= n )
        {
            if (it == 0){
                buffer[ep] = data[ep];
                bufferb[ep] = b[ep];
            }
            else
            {
                buffer[ep] = data[EPOCH*it-m+1+ep];
                bufferb[ep] = b[EPOCH*it-m+1+ep];
            }
            ep++;
        }

        /// Data are read in chunk of size EPOCH.
        /// When there is nothing to read, the loop ends. Could it happen that for a few points we dont compute the dtw?
        if (ep<=m-1)
            done = true;
        else
        {   lower_upper_lemire(buffer, ep, r, l_buff, u_buff);

            /// Just for printing a dot for approximate a million point. Not much accurate.
            if (it%(1000000/(EPOCH-m+1))==0)
                fprintf(stderr,".");

            /// Do main task here..
            ex=0; ex2=0;
            for(i=0; i<ep-1; i++) //ep-1 or ep?
            {

                /// A bunch of data has been read and pick one of them at a time to use
                d = buffer[i];

                /// Calcualte sum and sum square. need better accuracy maybe for floating point?
                ex += d;
                ex2 += d*d;
                curbit += bufferb[i];
                // printf(" %llu \n", i);
                // printf(" %3f \n", curbit);

                /// t is a circular array for keeping current data
                t[i%m] = d;
                /// Double the size for avoiding using modulo "%" operator
                t[(i%m)+m] = d;

                tb[i%m] = bufferb[i];
                /// Double the size for avoiding using modulo "%" operator
                tb[(i%m)+m] = bufferb[i];
                /// same for bitmask

                /// Start the task when there are more than m-1 points in the current chunk
                if( i >= m-1 )
                {
//                     fprintf(File3,"%d ",i);

                    // here an extra if: check if we're on a nan (bitmasked) pixel
                    mean = ex/m; //mean
                    std = ex2/m;
                    std = sqrt(std-mean*mean); //std

                    /// compute the start location of the data in the current circular array, t
                    j = (i+1)%m;
                    /// the start location of the data in the current chunk
                    I = i-(m-1);
                    // curbit=m;//tempdisable
                    if (curbit > m-0.5)
                    {
                    /// Use a constant lower bound to prune the obvious subsequence
                    lb_kim = lb_kim_hierarchy(t, q, j, m, mean, std, bsf);
                    if (lb_kim < bsf)
                    {
                        /// Use a linear time lower bound to prune; z_normalization of t will be computed on the fly.
                        /// uo, lo are envelop of the query.
                        lb_k = lb_keogh_cumulative(order, t, uo, lo, cb1, j, m, mean, std, bsf);
                        if (lb_k < bsf)
                        {
                            /// Take another linear time to compute z_normalization of t.
                            /// Note that for better optimization, this can merge to the previous function.
                            for(k=0;k<m;k++)
                            {   tz[k] = (t[(k+j)] - mean)/std;
                            }

                            /// Use another lb_keogh to prune
                            /// qo is the sorted query. tz is unsorted z_normalized data.
                            /// l_buff, u_buff are big envelop for all data in this chunk
                            lb_k2 = lb_keogh_data_cumulative(order, tz, qo, cb2, l_buff+I, u_buff+I, m, mean, std, bsf);
                            if (lb_k2 < bsf)
                            {
                                /// Choose better lower bound between lb_keogh and lb_keogh2 to be used in early abandoning DTW
                                /// Note that cb and cb2 will be cumulative summed here.
                                if (lb_k > lb_k2)
                                {
                                    cb[m-1]=cb1[m-1];
                                    for(k=m-2; k>=0; k--)
                                        cb[k] = cb[k+1]+cb1[k];
                                }
                                else
                                {
                                    cb[m-1]=cb2[m-1];
                                    for(k=m-2; k>=0; k--)
                                        cb[k] = cb[k+1]+cb2[k];
                                }

                                /// Compute DTW and early abandoning if possible
                                dist = dtw(tz, q, cb, m, r, bsf);

                                if( dist < bsf )
                                {   /// Update bsf
                                    /// loc is the real starting location of the nearest neighbor in the file
                                    bsf = dist;
                                    loc = (it)*(EPOCH-m+1) + i-m+1;
                                }
                            } else
                                keogh2++;
                        } else
                            keogh++;
                    } else
                        kim++;
                    } else
                        bm++;
                    /// Reduce obsolute points from sum and sum square
                    ex -= t[j];
                    ex2 -= t[j]*t[j];
                    curbit -= tb[j];
                }

            }

            /// If the size of last chunk is less then EPOCH, then no more data and terminate.
            if (ep<EPOCH)
                done=true;
            else
                it++;
            // done=true;
        }

        /// If the size of last chunk is less then EPOCH, then no more data and terminate.
        if (ep<EPOCH)
            done=true;
        else
            it++;


    }
//     fclose(File3);

    i= (it)*(EPOCH-m+1) + ep;

    *y = loc;
    *s = sqrt(bsf);

    t2 = clock();
    /// Note that loc and i are long long.
    // cout << "Location : " << loc << endl;
    // cout << "Distance : " << sqrt(bsf) << endl;
    // cout << "Data Scanned : " << i << endl;
    // cout << "Total Execution Time : " << (t2-t1)/CLOCKS_PER_SEC << " sec" << endl;
    // printf("\n");
    // printf("Pruned by bitmask    : %6.2f%%\n", ((double) bm / i)*100);
    // printf("Pruned by LB_Kim    : %6.2f%%\n", ((double) kim / i)*100);
    // printf("Pruned by LB_Keogh  : %6.2f%%\n", ((double) keogh / i)*100);
    // printf("Pruned by LB_Keogh2 : %6.2f%%\n", ((double) keogh2 / i)*100);
    // printf("DTW Calculation     : %6.2f%%\n", 100-(((double)kim+keogh+keogh2+bm)/i*100));
    // return

/*
    *y = loc;
    *s = sqrt(bsf);

    i= (it)*(EPOCH-m+1) + ep;

    free(q);
    free(u);
    free(l);
    free(uo);
    free(lo);
    free(qo);
    free(cb);
    free(cb1);
    free(cb2);
    free(tz);
    free(t);
    free(l_d);
    free(u_d);
    free(l_buff);
    free(u_buff);

    // return 0;


*/

}





//
// /// Main Function
// /// int main(  int argc , char *argv[] )
// void cpp_dtw(double *data, double *q, int n,int m, double R, double *y, double *s)
// {
//
//
//     double bsf;      /// best-so-far
//     double *t;       /// data array and query array
//     int *order;      ///new order of the query
//     double *u, *l, *qo, *uo, *lo,*tz,*cb, *cb1, *cb2,*u_d, *l_d;
//
//     // outputs..
//     *y = 0;
//     *s = 0;
//
//     /// name the parameters
//     double d;
//     long long i , j;
//     double ex , ex2 , mean, std;
//     int r = -1;
//     long long loc = 0;
//     double t1,t2;
//     int kim = 0,keogh = 0, keogh2 = 0;
//     double dist=0, lb_kim=0, lb_k=0, lb_k2=0;
//     double *buffer, *u_buff, *l_buff;
//     Index *Q_tmp;
//
//     /// For every EPOCH points, all cummulative values, such as ex (sum), ex2 (sum square), will be restarted for reducing the floating point error.
//     int EPOCH = 100000;
//
//     /// Sakoe-Chiba band
//     if (R<=1)
//             r = floor(R*m);
// 	else
//             r = floor(R);
//     // printf("%d\n",r);
//
//
//     /// start the clock
//     t1 = clock();
//
// // initializations
// {
//     qo = (double *)malloc(sizeof(double)*m);
//     uo = (double *)malloc(sizeof(double)*m);
//     lo = (double *)malloc(sizeof(double)*m);
//     order = (int *)malloc(sizeof(int)*m);
//     Q_tmp = (Index *)malloc(sizeof(Index)*m);
//     u = (double *)malloc(sizeof(double)*m);
//     l = (double *)malloc(sizeof(double)*m);
//     cb = (double *)malloc(sizeof(double)*m);
//     cb1 = (double *)malloc(sizeof(double)*m);
//     cb2 = (double *)malloc(sizeof(double)*m);
//     u_d = (double *)malloc(sizeof(double)*m);
//     l_d = (double *)malloc(sizeof(double)*m);
//     t = (double *)malloc(sizeof(double)*m*2);
//     tz = (double *)malloc(sizeof(double)*m);
//     buffer = (double *)malloc(sizeof(double)*EPOCH);
//     u_buff = (double *)malloc(sizeof(double)*EPOCH);
//     l_buff = (double *)malloc(sizeof(double)*EPOCH);
//
//     if( qo == NULL || uo == NULL || lo == NULL|| order == NULL|| Q_tmp == NULL || u == NULL || l == NULL || cb == NULL|| cb1 == NULL || cb2 == NULL ||u_d == NULL|| l_d == NULL || t == NULL ||  tz == NULL )
//     {
//         error(1);
//         return;
//     }
// }
//
//     /// Read query file
//     bsf = INF; /// So this finds global solution
//     j = 0;
//     ex = ex2 = 0;
//
//     // compute statistics for query /// alternatively use Error-free transformation here
//     for( i = 0 ; i < m ; i++ )
//     {
//         ex += q[i];    /// if the numbers are very big, would have some numerical issues?
//         ex2 += q[i]*q[i];
//     }
//     //Sum(a,b)
//
//     /// Do z-normalize the query, keep in same array, q
//     mean = ex/m;
//     std = ex2/m;
//     std = sqrt(std-mean*mean);
//     for( i = 0 ; i < m ; i++ )
//          q[i] = (q[i] - mean)/std;
//
//     lower_upper_lemire(q, m, r, l, u);
//
//     // printf("%f\n",l[0]);
//     // printf("%f\n",u[0]);
//     FILE *File2;
//     File2 = fopen("lemire_l.txt", "w+");
//     for(int k = 0 ; k < m ; k++ )
//         fprintf(File2,"%f ",l[k]);
//     fclose(File2);
//
//     File2 = fopen("lemire_u.txt", "w+");
//     for(int k = 0 ; k < m ; k++ )
//         fprintf(File2,"%f ",u[k]);
//     fclose(File2);
//
//
//     // test purposes: save lower and upper lemire to txt files
//
//     /// Sort the query one time by abs(z-norm(q[i]))
//     for( i = 0; i<m; i++)
//     {
//         Q_tmp[i].value = q[i];
//         Q_tmp[i].index = i;
//     }
//     qsort(Q_tmp, m, sizeof(Index),comp);
//
//     /// also create another arrays for keeping sorted envelop
//     for( i=0; i<m; i++)
//     {   int o = Q_tmp[i].index;
//         order[i] = o;
//         qo[i] = q[o];
//         uo[i] = u[o];
//         lo[i] = l[o];
//     }
//     free(Q_tmp);
//
//
//     /// Initial the cummulative lower bound / also consider uper bound maybe?
//     for( i=0; i<m; i++)
//     {   cb[i]=0;
//         cb1[i]=0;
//         cb2[i]=0;
//     }
//
//     i = 0;          /// current index of the data in current chunk of size EPOCH
//     j = 0;          /// the starting index of the data in the circular array, t
//     ex = ex2 = 0;
//     bool done = false;
//     int it=0, ep=0, k=0;
//     long long I;    /// the starting index of the data in current chunk of size EPOCH
//
//     /// Read first m-1 points
//     ep=0;
//     if (it==0)
//     {   for(k=0; k<m-1; k++)
//             buffer[k] = data[k];
//     }
//     printf("%f\n",buffer[0]);
//
// /*
//     while(!done)
//     {
//
//         /// Read first m-1 points
//         ep=0;
//         if (it==0)
//         {   for(k=0; k<m-1; k++)
//                 buffer[k] = data[k];
//         }
//         else
//         {   for(k=0; k<m-1; k++)
//                 buffer[k] = buffer[EPOCH-m+1+k];
//         }
//
//         /// Read buffer of size EPOCH or when all data has been read.
//         ep=m-1; // create buffer.. since all data is already loaded, might not be necessary?
//         while(ep<EPOCH && EPOCH*it+ep <= n )
//         {
//             buffer[ep] = data[EPOCH*it+ep];
//             ep++;
//         }
//
//         /// Data are read in chunk of size EPOCH.
//         /// When there is nothing to read, the loop is end. Could it happen that for a few points we dont compute the dtw?
//         if (ep<=m-1)
//         {   done = true;
//         } else
//         {   lower_upper_lemire(buffer, ep, r, l_buff, u_buff);
//
//             /// Just for printing a dot for approximate a million point. Not much accurate.
//             if (it%(1000000/(EPOCH-m+1))==0)
//                 fprintf(stderr,".");
//
//             /// Do main task here..
//             ex=0; ex2=0;
//             for(i=0; i<ep; i++)
//             {
//                 /// A bunch of data has been read and pick one of them at a time to use
//                 d = buffer[i];
//
//                 /// Calcualte sum and sum square. need better accuracy maybe for floating point?
//                 ex += d;
//                 ex2 += d*d;
//
//                 /// t is a circular array for keeping current data
//                 t[i%m] = d;
//
//                 /// Double the size for avoiding using modulo "%" operator
//                 t[(i%m)+m] = d;
//
//                 /// Start the task when there are more than m-1 points in the current chunk
//                 if( i >= m-1 )
//                 {
//                     mean = ex/m;
//                     std = ex2/m;
//                     std = sqrt(std-mean*mean);
//
//                     /// compute the start location of the data in the current circular array, t
//                     j = (i+1)%m;
//                     /// the start location of the data in the current chunk
//                     I = i-(m-1);
//
//                     /// Use a constant lower bound to prune the obvious subsequence
//                     lb_kim = lb_kim_hierarchy(t, q, j, m, mean, std, bsf);
//
//                     if (lb_kim < bsf)
//                     {
//                         /// Use a linear time lower bound to prune; z_normalization of t will be computed on the fly.
//                         /// uo, lo are envelop of the query.
//                         lb_k = lb_keogh_cumulative(order, t, uo, lo, cb1, j, m, mean, std, bsf);
//                         if (lb_k < bsf)
//                         {
//                             /// Take another linear time to compute z_normalization of t.
//                             /// Note that for better optimization, this can merge to the previous function.
//                             for(k=0;k<m;k++)
//                             {   tz[k] = (t[(k+j)] - mean)/std;
//                             }
//
//                             /// Use another lb_keogh to prune
//                             /// qo is the sorted query. tz is unsorted z_normalized data.
//                             /// l_buff, u_buff are big envelop for all data in this chunk
//                             lb_k2 = lb_keogh_data_cumulative(order, tz, qo, cb2, l_buff+I, u_buff+I, m, mean, std, bsf);
//                             if (lb_k2 < bsf)
//                             {
//                                 /// Choose better lower bound between lb_keogh and lb_keogh2 to be used in early abandoning DTW
//                                 /// Note that cb and cb2 will be cumulative summed here.
//                                 if (lb_k > lb_k2)
//                                 {
//                                     cb[m-1]=cb1[m-1];
//                                     for(k=m-2; k>=0; k--)
//                                         cb[k] = cb[k+1]+cb1[k];
//                                 }
//                                 else
//                                 {
//                                     cb[m-1]=cb2[m-1];
//                                     for(k=m-2; k>=0; k--)
//                                         cb[k] = cb[k+1]+cb2[k];
//                                 }
//
//                                 /// Compute DTW and early abandoning if possible
//                                 dist = dtw(tz, q, cb, m, r, bsf);
//
//                                 if( dist < bsf )
//                                 {   /// Update bsf
//                                     /// loc is the real starting location of the nearest neighbor in the file
//                                     bsf = dist;
//                                     loc = (it)*(EPOCH-m+1) + i-m+1;
//                                 }
//                             } else
//                                 keogh2++;
//                         } else
//                             keogh++;
//                     } else
//                         kim++;
//
//                     /// Reduce obsolute points from sum and sum square
//                     ex -= t[j];
//                     ex2 -= t[j]*t[j];
//                 }
//             }
//
//             /// If the size of last chunk is less then EPOCH, then no more data and terminate.
//             if (ep<EPOCH)
//                 done=true;
//             else
//                 it++;
//         }
//
//     }
//
//
//     i= (it)*(EPOCH-m+1) + ep;
//
//     free(q);
//     free(u);
//     free(l);
//     free(uo);
//     free(lo);
//     free(qo);
//     free(cb);
//     free(cb1);
//     free(cb2);
//     free(tz);
//     free(t);
//     free(l_d);
//     free(u_d);
//     free(l_buff);
//     free(u_buff);
//
//     t2 = clock();
// */
//
// }
//
//
//
//
// /*




   /// printf("\n");
/*
    /// Note that loc and i are long long.
    cout << "Location : " << loc << endl;
    cout << "Distance : " << sqrt(bsf) << endl;
    cout << "Data Scanned : " << i << endl;
    cout << "Total Execution Time : " << (t2-t1)/CLOCKS_PER_SEC << " sec" << endl;

    */
    /*
    *y = loc;
    *s = sqrt(bsf);

    /// printf is just easier for formating ;)
    ///*
    printf("\n");
    printf("Pruned by LB_Kim    : %6.2f%%\n", ((double) kim / i)*100);
    printf("Pruned by LB_Keogh  : %6.2f%%\n", ((double) keogh / i)*100);
    printf("Pruned by LB_Keogh2 : %6.2f%%\n", ((double) keogh2 / i)*100);
    printf("DTW Calculation     : %6.2f%%\n", 100-(((double)kim+keogh+keogh2)/i*100));
   //*/
    /// return 0;

// }



/* The gateway function that replaces the "main".
 *plhs[] - the array of output values
 *prhs[] - the array of input values*/

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    /// output location
    double *y;
    /// output score
    double *s;

	double R;
    /*Variable declarations as in C*/
    ///  char Data_File, Query_File;
    int lenQ;
    int lenD;


    ///int L; /// length of the overlap

	int buflen0,status0;
    int buflen1,status1;

	char *input_buf0, *output_buf0;
	char *input_buf1, *output_buf1;

 ///
 ///   Query_File = mxGetString(prhs[1]);


    /* Check for proper number of arguments. */
    if (nrhs != 6) {
        mexErrMsgTxt("Six inputs required.");
	}
    else if (nlhs > 2) {
        mexErrMsgTxt("Too many output arguments");
    }

    lenD = mxGetScalar(prhs[3]); /// data length
    lenQ = mxGetScalar(prhs[4]); /// query length
    R = mxGetScalar(prhs[5]); /// Sakoe-Chiba width

    //   /* Get the length of the input string. */
    // buflen0 = (mxGetM(prhs[0]) * mxGetN(prhs[0])) + 1;
    // /* Allocate memory for input and output strings. */
	// input_buf0 =(char *) mxCalloc(buflen0, sizeof(char));
    //
    // /* Copy the string data from prhs[0] into a C string
    // * input_buf. */
   	// status0 = mxGetString(prhs[0], input_buf0, buflen0);

	// /* Get the length of the input string. */
    // buflen1 = (mxGetM(prhs[1]) * mxGetN(prhs[1])) + 1;
    // /* Allocate memory for input and output strings. */
	// input_buf1 =(char *) mxCalloc(buflen1, sizeof(char));

    /* Copy the string data from prhs[0] into a C string
    * input_buf. */
   	// status1 = mxGetString(prhs[1], input_buf1, buflen1);

    // here we load both data and query..
    double *d;
    d = (double *)malloc(sizeof(double)*lenD);
    d = mxGetPr(prhs[0]);

    double *q;
    q = (double *)malloc(sizeof(double)*lenQ);
    q = mxGetPr(prhs[1]);

    double *b;
    b = (double *)malloc(sizeof(double)*lenD);
    b = mxGetPr(prhs[2]);
    // b = (int *) mxGetPr(prhs[2]);

    // print a single value
    // mexPrintf("%4.3f\n", b[99]);
    // mexPrintf("%4.3f\n", b[100]);
    // mexPrintf("%4.3f\n", b[101]);

    // mexPrintf("%f\n", d[0]);

    // q =  mxGetData(prhs[1]);
    // mexPrintf("%f\n", q[0]);
//     /*
    // mexPrintf("Number of inputs:  %d\n", nrhs);
    // mexPrintf("Number of outputs: %d\n", nlhs);
	// mexPrintf("lenQ: %d\n", lenQ);
    // mexPrintf("lenD: %d\n", lenD);
    // mexPrintf("R: %4.3f\n", R);
//     */

	plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
	y = mxGetPr(plhs[0]);

    plhs[1] = mxCreateDoubleMatrix(1, 1, mxREAL);
	s = mxGetPr(plhs[1]);


    cpp_dtw(d, q,b, lenD, lenQ, R, y, s);


}
