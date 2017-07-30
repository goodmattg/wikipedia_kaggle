/***********************************************************************/
/************************* DISCLAIMER **********************************/
/***********************************************************************/
/** This UCR Suite software is copyright protected ? 2012 by          **/
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

/** Heavily modified by Matthew Goodman for nearest neighbor          **/
/**  among set of sequences                                           **/

/***********************************************************************/
/***********************************************************************/



#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <iostream>
#include <cmath>

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
void lower_upper_lemire(double *t, int len, int offset, int r, double *l, double *u)
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
            u[i-r-1] = t[offset + front(&du)];
            l[i-r-1] = t[offset + front(&dl)];
        }
        if (t[offset + i] > t[offset + i-1])
        {
            pop_back(&du);
            while (!empty(&du) && t[offset + i] > t[offset + back(&du)])
                pop_back(&du);
        }
        else
        {
            pop_back(&dl);
            while (!empty(&dl) && t[offset + i] < t[offset + back(&dl)])
                pop_back(&dl);
        }
        push_back(&du, i);
        push_back(&dl, i);
        if (i == 2 * r + 1 + front(&du))
            pop_front(&du);
        else if (i == 2 * r + 1 + front(&dl))
            pop_front(&dl);
    }
    // POSSIBLE BUG?? len+r-1   ??
    for (int i = len; i < len+r+1; i++)
    {
        u[i-r-1] = t[offset + front(&du)];
        l[i-r-1] = t[offset + front(&dl)];
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
double lb_kim_hierarchy(double *t, double *q, int j, int len, double bsf = INF)
{

    // ALREADY NORMALIZED

    /// 1 point at front and back
    double d, lb;
    double x0 = t[j];
    double y0 = t[(len-1+j)];
    lb = dist(x0,q[0]) + dist(y0,q[len-1]);
    if (lb >= bsf)   return lb;

    /// 2 points at front
    double x1 = t[(j+1)];
    d = min(dist(x1,q[0]), dist(x0,q[1]));
    d = min(d, dist(x1,q[1]));
    lb += d;
    if (lb >= bsf)   return lb;

    /// 2 points at back
    double y1 = t[(len-2+j)];
    d = min(dist(y1,q[len-1]), dist(y0, q[len-2]) );
    d = min(d, dist(y1,q[len-2]));
    lb += d;
    if (lb >= bsf)   return lb;

    /// 3 points at front
    double x2 = t[(j+2)];
    d = min(dist(x0,q[2]), dist(x1, q[2]));
    d = min(d, dist(x2,q[2]));
    d = min(d, dist(x2,q[1]));
    d = min(d, dist(x2,q[0]));
    lb += d;
    if (lb >= bsf)   return lb;

    /// 3 points at back
    double y2 = t[(len-3+j)];
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
/// u_query_o, l_query_o: upper and lower envelops for the query, which already sorted.
/// t     : array keeping the current data.
/// j     : index of the starting location in t
/// cb    : (output) current bound at each position. It will be used later for early abandoning in DTW.
double lb_keogh_cumulative(int* order, double *t, double *u_query_o, double *l_query_o,
                           double *cb, int j, int len,
                           double best_so_far = INF)
{
    double lb = 0;
    double x, d;

    for (int i = 0; i < len && lb < best_so_far; i++)
    {
        // Data already mean normalized
        x = t[(order[i]+j)];
        d = 0;
        if (x > u_query_o[i])
            d = dist(x, u_query_o[i]);
        else if(x < l_query_o[i])
            d = dist(x, l_query_o[i]);
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
double lb_keogh_data_cumulative(int* order, double *qo, double *cb,
                                double *l_data, double *u_data, int len,
                                double best_so_far = INF)
{
    double lb = 0;
    double uu,ll,d;

    for (int i = 0; i < len && lb < best_so_far; i++)
    {
        uu = u_data[order[i]];
        ll = l_data[order[i]];

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
double dtw(double* A, int offset, double* B, double *cb, int m, int r, double bsf = INF)
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
                cost[k]=dist(A[offset+0],B[0]);
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
            cost[k] = min( min( x, y) , z) + dist(A[offset+i],B[j]);

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
        printf("Command Usage:  UCR_DTW.exe  data-file  query-file   m   R output-file\n\n");
        printf("For example  :  UCR_DTW.exe  data.txt   query.txt   128  0.05 output.txt\n");
    }
    else if (id == 5)
    {
        printf("LB_Kim reeturn NaN");
    }

    fclose(query_fptr);
    fclose(data_fptr);
    fclose(out_fptr);

    exit(1);
}




/// ------------------------------------------------------------------------------

/// Main Function
// int main(  int argc , char *argv[] )
// {
//     FILE *fp;            /// data file pointer
//     FILE *qp;            /// query file pointer
//     double bsf;          /// best-so-far
//     double *q;       /// data array and query array
//     int *order;          ///new order of the query
//     double *qo, *u_query, *l_query,  *u_query_o, *l_query_o,*cb, *cb1, *cb2;

//     double d;
//     long long i , j;
//     double ex , ex2 , mean, std;
//     int m=-1, r=-1;
//     long long loc = 0;
//     double t1,t2;
//     int kim = 0,keogh = 0, keogh2 = 0;
//     double dist=0, lb_kim=0, lb_k=0, lb_k2=0;
//     double *buffer, *u_buff, *l_buff;
//     Index *Q_tmp;


//     int EPOCH;
//     int chunks_in_epoch; // "chunks in epoch"
//     double *epoch_ex, *epoch_ex2, *epoch_mean, *epoch_std;

//     /// If not enough input, display an error.
//     if (argc<=3)
//         error(4);

//     /// read size of the query
//     if (argc>3)
//         m = atol(argv[3]);

//     /// For every EPOCH points, all cummulative values, such as ex (sum), ex2 (sum square), will be restarted for reducing the floating point error.
//     /// EPOCH set as multiple of m
//     if (m > 100000) {
//         chunks_in_epoch = 1;
//         EPOCH = m;
//     } else {
//         chunks_in_epoch = 100000/m;
//         EPOCH = chunks_in_epoch * m;
//     }
//     printf("Chunks Per Epoch: %d\n",chunks_in_epoch);

//     /// read warping windows
//     if (argc>4)
//     {   double R = atof(argv[4]);
//         if (R<=1)
//             r = floor(R*m);
//         else
//             r = floor(R);
//     }
//     printf("Warping window width: %d\n", r);

//     fp = fopen(argv[1],"r");
//     if( fp == NULL )
//         error(2);

//     qp = fopen(argv[2],"r");
//     if( qp == NULL )
//         error(2);

//     /// start the clock
//     t1 = clock();

//     /**
//     ==================================================================
//     **/

//     // Query Data
//     q = (double *)malloc(sizeof(double)*m);
//     if( q == NULL )
//         error(1);

//     // Ordered Query_Data
//     qo = (double *)malloc(sizeof(double)*m);
//     if( qo == NULL )
//         error(1);

//     // Query Upper Envelope
//     u_query = (double *)malloc(sizeof(double)*m);
//     if( u_query == NULL )
//         error(1);

//     // Ordered Query Upper Envelope
//     u_query_o = (double *)malloc(sizeof(double)*m);
//     if( u_query_o == NULL )
//         error(1);

//     // Query Lower Envelope
//     l_query = (double *)malloc(sizeof(double)*m);
//     if( l_query == NULL )
//         error(1);

//     //  Ordered Query Lower Envelope
//     l_query_o = (double *)malloc(sizeof(double)*m);
//     if( l_query_o == NULL )
//         error(1);

//     // Indices of sorted query data
//     order = (int *)malloc(sizeof(int)*m);
//     if( order == NULL )
//         error(1);

//     // Temporary Structure for query sorting
//     Q_tmp = (Index *)malloc(sizeof(Index)*m);
//     if( Q_tmp == NULL )
//         error(1);

//     buffer = (double *)malloc(sizeof(double)*EPOCH);
//     if( buffer == NULL )
//         error(1);

//     ///------------------------------------------------------------

//     epoch_ex = (double *)malloc(sizeof(double)*chunks_in_epoch);
//     if( epoch_ex == NULL )
//         error(1);

//     epoch_ex2 = (double *)malloc(sizeof(double)*chunks_in_epoch);
//     if( epoch_ex2 == NULL )
//         error(1);

//     epoch_mean = (double *)malloc(sizeof(double)*chunks_in_epoch);
//     if( epoch_mean == NULL )
//         error(1);

//     epoch_std = (double *)malloc(sizeof(double)*chunks_in_epoch);
//     if( epoch_std == NULL )
//         error(1);

//     u_buff = (double *)malloc(sizeof(double)*m);
//     if( u_buff == NULL )
//         error(1);

//     l_buff = (double *)malloc(sizeof(double)*m);
//     if( l_buff == NULL )
//         error(1);

//     cb = (double *)malloc(sizeof(double)*m);
//     if( cb == NULL )
//         error(1);

//     cb1 = (double *)malloc(sizeof(double)*m);
//     if( cb1 == NULL )
//         error(1);

//     cb2 = (double *)malloc(sizeof(double)*m);
//     if( cb2 == NULL )
//         error(1);


//     /// Read query file
//     bsf = INF;
//     i = 0;
//     j = 0;
//     ex = ex2 = 0;

//     while(fscanf(qp,"%lf",&d) != EOF && i < m)
//     {
//         ex += d;
//         ex2 += d*d;
//         q[i] = d; // q is the query file data
//         i++;
//     }
//     fclose(qp);

//     /// Z-normalize the query inplace
//     mean = ex/m;
//     std = ex2/m;
//     std = sqrt(std-mean*mean);
//     for( i = 0 ; i < m ; i++ ) {
//         q[i] = (q[i] - mean)/std;
//         if (isnan(q[i])) {
//             q[i] = 0.0;
//         }
//     }


//     /// Create envelop of the query
//     /// lower envelop:  l_query
//     /// upper envelop:  u_query
//     lower_upper_lemire(q, m, 0, r, l_query, u_query);

//     /// Sort the query one time by abs(z-norm(q[i]))
//     for( i = 0; i<m; i++)
//     {
//         Q_tmp[i].value = q[i];
//         Q_tmp[i].index = i;
//     }
//     qsort(Q_tmp, m, sizeof(Index),comp);
//     cout << "Query sorted\n";

//     /// also create another arrays for keeping sorted envelop
//     for(i=0; i<m; i++)
//     {   int o = Q_tmp[i].index;
//         order[i] = o;
//         qo[i] = q[o];
//         u_query_o[i] = u_query[o];
//         l_query_o[i] = l_query[o];
//     }
//     free(Q_tmp);

//     /// Initialize the cummulative lower bound arrays
//     for( i=0; i<m; i++)
//     {   cb[i]=0;
//         cb1[i]=0;
//         cb2[i]=0;
//     }


//     i = 0;  /// current index of the starting pt. in EPOCH -> [i, i+m-1]
//     ex = ex2 = 0;
//     bool done = false;
//     int it=0, k=0;
//     int vals_in_epoch;
//     int epoch_count = 0;

//     while(!done)
//     {
//          printf("\nStarting Epoch : %d\n", it);

//         /// Attempt to read in EPOCH and accumulate data for later z-normalize
//         int chunk_ctr;
//         vals_in_epoch = 0;
//         for(int k=0; k<EPOCH; k++)
//         {
//             chunk_ctr = k/m;
//             if (fscanf(fp,"%lf",&d) != EOF)
//             {
//                 buffer[k] = d;
//                 epoch_ex[chunk_ctr] += d;
//                 epoch_ex2[chunk_ctr] += d*d;
//             } else
//             {
//                 printf("Hit EOF while reading data\n");
//                 if (k < m-1)
//                 {
//                     // error here
//                     cout << "weird values read in less than m-1";
//                     break;
//                 } else
//                 {
//                     vals_in_epoch = k;
//                     done = true;
//                     break;
//                 }

//             }
//         }
//         // Number of values in the epoch
//         if (!vals_in_epoch) {
//             vals_in_epoch = EPOCH;
//         }
//         printf("Vals in epoch: %d\n", vals_in_epoch);

//         // Z-normalize routine
//         for (int ch =0 ; ch < chunks_in_epoch ; ch++) {
//             // mean = ex/m;
//             // std = ex2/m;
//             // std = sqrt(std-mean*mean);
//             epoch_mean[ch] = epoch_ex[ch]/m;
//             epoch_std[ch] = epoch_ex2[ch]/m;
//             epoch_std[ch] = sqrt(epoch_std[ch]-epoch_mean[ch]*epoch_mean[ch]);
//             // Now Z-normalize each chunk in the buffer
//             for (int idx = ch*m ; idx < (ch+1)*m ; idx++) {
//                 double tmp = buffer[idx];
//                 buffer[idx] = (buffer[idx] - epoch_mean[ch]) / epoch_std[ch];
//                 if (isnan(buffer[idx])) {
//                     // Weird case where entire row is zero - just return 0.0 for all values
//                     buffer[idx] = 0.0;
//                 }
//             }
//         }
//         printf("Epoch z-normalized\n");

//         // Now iterate over chunks in epoch to compute DTW
//         for (int ch = 0 ; ch < chunks_in_epoch ; ch++) {

//             // Create lower-envelope and upper-envolope of the chunk
//             lower_upper_lemire(buffer, m, ch*m, r, l_buff, u_buff);
//             // l_buff   :   chunk lower envelope
//             // u_buff   :   chunk uppper envelope

//             /// Just for printing a dot for approximate a million point.
//             fprintf(stderr,".");
//             /// Use a constant lower bound to prune the obvious subsequence. VERIFIED

//             if (isnan(lb_kim = lb_kim_hierarchy(buffer, q, ch*m, m, bsf)))
//                 error(5);

//             if (lb_kim < bsf)
//             {
//                 // Use a linear time lower bound to prune;
//                 // u_query, lo are envelop of the query.
//                 // VERIFIED
//                 lb_k = lb_keogh_cumulative(order, buffer, u_query_o, l_query_o, cb1, ch*m, m, bsf);

//                 if (lb_k < bsf)
//                 {
//                     /// Use another lb_keogh to prune
//                     /// qo is the sorted query. tz is unsorted z_normalized data.
//                     /// l_buff, u_buff are big envelop for all data in this chunk
//                     lb_k2 = lb_keogh_data_cumulative(order, qo, cb2, l_buff, u_buff, m, bsf);
//                     if (lb_k2 < bsf)
//                     {
//                         /// Choose better lower bound between lb_keogh and lb_keogh2 to be used in early abandoning DTW
//                         /// Note that cb and cb2 will be cumulative summed here.
//                         if (lb_k > lb_k2)
//                         {
//                             cb[m-1]=cb1[m-1];
//                             for(k=m-2; k>=0; k--)
//                                 cb[k] = cb[k+1]+cb1[k];
//                         }
//                         else
//                         {
//                             cb[m-1]=cb2[m-1];
//                             for(k=m-2; k>=0; k--)
//                                 cb[k] = cb[k+1]+cb2[k];
//                         }

//                         /// Compute DTW and early abandoning if possible
//                         dist = dtw(buffer, ch*m, q, cb, m, r, bsf);

//                         if( dist < bsf )
//                         {   /// Update bsf
//                             /// loc is the real starting location of the nearest neighbor in the file
//                             bsf = dist;
//                             loc = it;
//                             printf("BSF val: %f\n",bsf);
//                         }
//                     } else
//                         keogh2++;
//                 } else
//                     keogh++;
//             } else
//                 kim++;

//             it++;
//         }
//         /// If the size of last chunk is less then EPOCH, then no more data and terminate.
//         epoch_count++;
//     }


//     // i = (it)*(EPOCH-m+1) + ep;
//     fclose(fp);
//     free(q);
//     free(u_query);
//     free(l_query);
//     free(l_query_o);
//     free(u_query_o);
//     free(cb);
//     free(cb1);
//     free(cb2);
//     free(l_buff);
//     free(u_buff);

//     t2 = clock();
//     printf("\n");

//     printf("Location: %lld (row %lld)\n", loc, loc+1);
//     cout << "Distance : " << sqrt(bsf) << endl;
//     cout << "Total Execution Time : " << (t2-t1)/CLOCKS_PER_SEC << " sec" << endl;

//     /// printf is just easier for formating ;)
//     printf("\n");
//     printf("Pruned by LB_Kim    : %6.2f%%\n", ((double) kim / it)*100);
//     printf("Pruned by LB_Keogh  : %6.2f%%\n", ((double) keogh / it)*100);
//     printf("Pruned by LB_Keogh2 : %6.2f%%\n", ((double) keogh2 / it)*100);
//     printf("DTW Calculation     : %6.2f%%\n", 100-(((double)kim+keogh+keogh2)/it*100));
//     return 0;
// }










































int NearestNeighborSearch(FILE *data_fptr, FILE *query_fptr, FILE *out_fptr, int m, int r) {

    double bsf;          /// best-so-far
    double *q;       /// data array and query array
    int *order;          ///new order of the query
    double *qo, *u_query, *l_query,  *u_query_o, *l_query_o,*cb, *cb1, *cb2;

    double d;
    long long i , j;
    double ex , ex2 , mean, std;
    long long loc = 0;
    double t1, t2;
    int kim = 0,keogh = 0, keogh2 = 0;
    double dist=0, lb_kim=0, lb_k=0, lb_k2=0;
    double *buffer, *u_buff, *l_buff;
    Index *Q_tmp;

    // RUN-TIME FLAGS
    bool QUERY_FINISHED = false;
    bool DATA_FINISHED = false;

    int EPOCH;
    int chunks_in_epoch; // "chunks in epoch"
    double *epoch_ex, *epoch_ex2, *epoch_mean, *epoch_std;

    /*
    For every EPOCH points, all cummulative values, such as ex (sum),
    ex2 (sum square), will be restarted for reducing the floating point error.
    EPOCH set as multiple of m
    */
    if (m > 100000) {
        chunks_in_epoch = 1;
        EPOCH = m;
    } else {
        chunks_in_epoch = 100000/m;
        EPOCH = chunks_in_epoch * m;
    }
    printf("Chunks Per Epoch: %d\n",chunks_in_epoch);


    t1 = clock();

    // Query Data
    q = (double *)malloc(sizeof(double)*m);
    if( q == NULL )
        error(1);

    // Ordered Query_Data
    qo = (double *)malloc(sizeof(double)*m);
    if( qo == NULL )
        error(1);

    // Query Upper Envelope
    u_query = (double *)malloc(sizeof(double)*m);
    if( u_query == NULL )
        error(1);

    // Ordered Query Upper Envelope
    u_query_o = (double *)malloc(sizeof(double)*m);
    if( u_query_o == NULL )
        error(1);

    // Query Lower Envelope
    l_query = (double *)malloc(sizeof(double)*m);
    if( l_query == NULL )
        error(1);

    //  Ordered Query Lower Envelope
    l_query_o = (double *)malloc(sizeof(double)*m);
    if( l_query_o == NULL )
        error(1);

    // Indices of sorted query data
    order = (int *)malloc(sizeof(int)*m);
    if( order == NULL )
        error(1);

    // Temporary Structure for query sorting
    Q_tmp = (Index *)malloc(sizeof(Index)*m);
    if( Q_tmp == NULL )
        error(1);

    buffer = (double *)malloc(sizeof(double)*EPOCH);
    if( buffer == NULL )
        error(1);

    epoch_ex = (double *)malloc(sizeof(double)*chunks_in_epoch);
    if( epoch_ex == NULL )
        error(1);

    epoch_ex2 = (double *)malloc(sizeof(double)*chunks_in_epoch);
    if( epoch_ex2 == NULL )
        error(1);

    epoch_mean = (double *)malloc(sizeof(double)*chunks_in_epoch);
    if( epoch_mean == NULL )
        error(1);

    epoch_std = (double *)malloc(sizeof(double)*chunks_in_epoch);
    if( epoch_std == NULL )
        error(1);

    u_buff = (double *)malloc(sizeof(double)*m);
    if( u_buff == NULL )
        error(1);

    l_buff = (double *)malloc(sizeof(double)*m);
    if( l_buff == NULL )
        error(1);

    cb = (double *)malloc(sizeof(double)*m);
    if( cb == NULL )
        error(1);

    cb1 = (double *)malloc(sizeof(double)*m);
    if( cb1 == NULL )
        error(1);

    cb2 = (double *)malloc(sizeof(double)*m);
    if( cb2 == NULL )
        error(1);

    while (!QUERY_FINISHED)

    {
        /// Line in the query file
        bsf = INF;
        i = j = 0;
        ex = ex2 = 0;

        for (i = 0; i < m; i++) {
            if (fscanf(query_fptr,"%lf",&d) != EOF) {
                ex += d;
                ex2 += d*d;
                q[i] = d;
            } else {
                QUERY_FINISHED = true;
            }
        }
        fprintf(stderr,".\n");


        /// Z-normalize the query inplace
        mean = ex/m;
        std = ex2/m;
        std = sqrt(std-mean*mean);
        for( i = 0 ; i < m ; i++ ) {
            q[i] = (q[i] - mean)/std;
            if (isnan(q[i])) {
                q[i] = 0.0;
            }
        }

        /// Create envelop of the query
        /// lower envelop:  l_query
        /// upper envelop:  u_query
        lower_upper_lemire(q, m, 0, r, l_query, u_query);

        /// Sort the query one time by abs(z-norm(q[i]))
        for( i = 0; i<m; i++)
        {
            Q_tmp[i].value = q[i];
            Q_tmp[i].index = i;
        }
        qsort(Q_tmp, m, sizeof(Index),comp);
        // cout << "Query sorted\n";

        /// also create another arrays for keeping sorted envelop
        for(i=0; i<m; i++)
        {   int o = Q_tmp[i].index;
            order[i] = o;
            qo[i] = q[o];
            u_query_o[i] = u_query[o];
            l_query_o[i] = l_query[o];
        }

        /// Initialize the cummulative lower bound arrays
        for( i=0; i<m; i++)
        {   cb[i]=0;
            cb1[i]=0;
            cb2[i]=0;
        }

        i = 0;  /// current index of the starting pt. in EPOCH -> [i, i+m-1]
        ex = ex2 = 0;
        int it=0, k=0;
        int vals_in_epoch;
        int epoch_count = 0;


        DATA_FINISHED = false;
        rewind(data_fptr);
        // printf("Starting on Data");


        while(!DATA_FINISHED)
        {
            // printf("\nStarting Epoch : %d\n", it);
            /// Attempt to read in EPOCH and accumulate data for later z-normalize
            int chunk_ctr;
            vals_in_epoch = 0;
            for(int k=0; k<EPOCH; k++)
            {
                chunk_ctr = k/m;
                if (fscanf(data_fptr,"%lf",&d) != EOF)
                {
                    buffer[k] = d;
                    epoch_ex[chunk_ctr] += d;
                    epoch_ex2[chunk_ctr] += d*d;
                } else
                {
                    // printf("Hit EOF while reading data\n");
                    if (k < m-1)
                    {
                        // error here
                        cout << "weird values read in less than m-1";
                        break;
                    } else
                    {
                        vals_in_epoch = k;
                        DATA_FINISHED = true;
                        break;
                    }
                }
            }
            // Number of values in the epoch
            if (!vals_in_epoch) {
                vals_in_epoch = EPOCH;
            }
            // printf("Vals in epoch: %d\n", vals_in_epoch);

            // Z-normalize routine
            for (int ch =0 ; ch < chunks_in_epoch ; ch++) {
                epoch_mean[ch] = epoch_ex[ch]/m;
                epoch_std[ch] = epoch_ex2[ch]/m;
                epoch_std[ch] = sqrt(epoch_std[ch]-epoch_mean[ch]*epoch_mean[ch]);
                // Now Z-normalize each chunk in the buffer
                for (int idx = ch*m ; idx < (ch+1)*m ; idx++) {
                    buffer[idx] = (buffer[idx] - epoch_mean[ch]) / epoch_std[ch];
                    if (isnan(buffer[idx])) {
                        // Weird case where entire row is zero - just return 0.0 for all values
                        buffer[idx] = 0.0;
                    }
                }
            }
            // printf("Epoch z-normalized\n");

            // Now iterate over chunks in epoch to compute DTW
            for (int ch = 0 ; ch < chunks_in_epoch ; ch++) {

                // Create lower-envelope and upper-envolope of the chunk
                lower_upper_lemire(buffer, m, ch*m, r, l_buff, u_buff);
                // l_buff   :   chunk lower envelope
                // u_buff   :   chunk uppper envelope

                /// Just for printing a dot for approximate a million point.
                /// Use a constant lower bound to prune the obvious subsequence. VERIFIED

                if (isnan(lb_kim = lb_kim_hierarchy(buffer, q, ch*m, m, bsf)))
                    error(5);

                if (lb_kim < bsf)
                {
                    // Use a linear time lower bound to prune;
                    // u_query, lo are envelop of the query.
                    // VERIFIED
                    lb_k = lb_keogh_cumulative(order, buffer, u_query_o, l_query_o, cb1, ch*m, m, bsf);

                    if (lb_k < bsf)
                    {
                        /// Use another lb_keogh to prune
                        /// qo is the sorted query. tz is unsorted z_normalized data.
                        /// l_buff, u_buff are big envelop for all data in this chunk
                        lb_k2 = lb_keogh_data_cumulative(order, qo, cb2, l_buff, u_buff, m, bsf);
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
                            dist = dtw(buffer, ch*m, q, cb, m, r, bsf);

                            if( dist < bsf )
                            {   /// Update bsf
                                /// loc is the real starting location of the nearest neighbor in the file
                                bsf = dist;
                                loc = it;
                                // printf("BSF val: %f\n",bsf);
                            }
                        } else
                            keogh2++;
                    } else
                        keogh++;
                } else
                    kim++;

                it++;
            }
            /// If the size of last chunk is less then EPOCH, then no more data and terminate.
            epoch_count++;
        }

        // printf("Data Round done\n");
        fprintf(out_fptr, "%lld\n", loc);

        // printf("Location: %lld (row %lld)\n", loc, loc+1);
        // cout << "Distance : " << sqrt(bsf) << endl;

    }

        fclose(query_fptr);
        fclose(data_fptr);
        fclose(out_fptr);



        free(Q_tmp);
        free(q);
        free(u_query);
        free(l_query);
        free(l_query_o);
        free(u_query_o);
        free(l_buff);
        free(u_buff);
        free(cb);
        free(cb1);
        free(cb2);

        t2 = clock();
        cout << "Total Execution Time : " << (t2-t1)/CLOCKS_PER_SEC << " sec" << endl;
        printf("\n");
        return 0;

}
























int main(int argc, char *argv[]) {
    // Data_File            | argv[1]
    // Query_File           | argv[2]
    // M                    | argv[3]
    // R                    | argv[4]
    // Search_Output_File   | argv[5]

    FILE *d_fptr; // data file pointer
    FILE *q_fptr; // query file pointer
    FILE *o_fptr;  // output file pointer
    int m=-1, r=-1;

    if (argc < 5) {
        error(4);
    }

    // read sequence length
    m = atol(argv[3]);

    // read warping windows
    double R = atof(argv[4]);
    r = (R<=1) ? floor(R*m) : floor(R);

    printf("Warping window width: %d\n", r);

    // Open data file for reading
    d_fptr = fopen(argv[1], "r");
    if( d_fptr == NULL )
        error(2);

    // Open query file for reading
    q_fptr = fopen(argv[2], "r");
    if( q_fptr == NULL )
        error(2);

    // Open output file for writing
    o_fptr = fopen(argv[5], "w");
    if( o_fptr == NULL )
        error(2);


    /* ---------------------------------------------------------------------- */
    // Run the nearest neighor search
    NearestNeighborSearch(d_fptr, q_fptr, o_fptr, m, r);

    /* ---------------------------------------------------------------------- */

    return 0;

}



