#include "fasp.h"
#include "fasp_functs.h"

#if 0
void generate_S_theta ( dCSRmat *A, 
                        iCSRmat *S, 
                        REAL theta )
{
    const INT row=A->row, col=A->col;
    const INT row_plus_one = row+1;
    const INT nnz=A->IA[row]-A->IA[0];
    
    INT index, i, j, begin_row, end_row;
    INT *ia=A->IA, *ja=A->JA;
    REAL *aj=A->val;

    // get the diagnal entry of A
    //dvector diag; fasp_dcsr_getdiag(0, A, &diag);
    
    /* generate S */
    REAL row_abs_sum;     

    // copy the structure of A to S
    S->row=row; S->col=col; S->nnz=nnz; S->val=NULL;
    
    S->IA=(INT*)fasp_mem_calloc(row_plus_one, sizeof(INT));
    
    S->JA=(INT*)fasp_mem_calloc(nnz, sizeof(INT));
    
    fasp_iarray_cp(row_plus_one, ia, S->IA);
    fasp_iarray_cp(nnz, ja, S->JA);
    
    for (i=0;i<row;++i) {
        /* compute scaling factor and row sum */
        row_abs_sum=0;
            
        begin_row=ia[i]; end_row=ia[i+1];

        for (j=begin_row;j<end_row;j++) {
            row_abs_sum+=ABS(aj[j]);
        }
	row_abs_sum = row_abs_sum*theta;

        /* deal with the diagonal element of S */
//        for (j=begin_row;j<end_row;j++) {
//            if (ja[j]==i) {S->JA[j]=-1; break;}
//        }
            
        /* deal with  the element of S */
        for (j=begin_row;j<end_row;j++){
        /* if $\sum_{j=1}^n |a_{ij}|*theta>= |a_{ij}|$ */
            //if ((row_abs_sum >= ABS(aj[j]))||(ja[j]==i)) {		    
            if ( (row_abs_sum >= ABS(aj[j])) && (ja[j] !=i) )  {		    
                S->JA[j]=-1;
//		printf("S->JA[%d]=%d,row_abs_sum = %f, ABS(aj[j]) = %f\n",j,S->JA[j],row_abs_sum,ABS(aj[j]));
	    } 
        }
    } // end for i
    
    /* Compress the strength matrix */
    index=0;
    for (i=0;i<row;++i) {
        S->IA[i]=index;
        begin_row=ia[i]; end_row=ia[i+1]-1;
        for (j=begin_row;j<=end_row;j++) {
            if (S->JA[j]>-1) {
                S->JA[index]=S->JA[j];
                index++;
//		printf("S->JA[%d]=%d\n",index,S->JA[index]);
            }
        }
    }
    
    if (index > 0) {
        S->IA[row]=index;
        S->nnz=index;
        S->JA=(INT*)fasp_mem_realloc(S->JA,index*sizeof(INT));
    }
    else {
        S->nnz = 0;
        S->JA = NULL;
    }
    
}


void dCSRmat_Multicoloring_Theta(AMG_data *mgl,
				                 REAL theta,
                                 INT *rowmax,
                                 INT *groups)
{
#if MULTI_COLORS_ORDER    
    INT k,i,j,pre,group;
    INT igold,iend,iavg;
    INT icount;
    INT front,rear;
    dCSRmat A = mgl->A;
    INT n = A.row;
     iCSRmat S;

     INT *IA,*JA;

     if (theta > 0 && theta < 1.0) {
     generate_S_theta(&A, &S, theta);
     IA = S.IA;
     JA = S.JA;
     } else if (theta == 1.0 ){

     mgl->ic = (INT *)malloc(sizeof(INT)*2);
     mgl->icmap = (INT *)malloc(sizeof(INT)*(n+1));
     mgl->ic[0] = 0;
     mgl->ic[1] = n;
     for(k=0; k<n; k++)  mgl->icmap[k]= k;

     mgl->colors = 1;
     *groups = 1;
     *rowmax = 1;

     printf( "Theta = %lf \n",theta );

     return;

     } else{
     IA = A.IA;
     JA = A.JA;
     }
//---------------------------------------------------------------------------
    INT *cq = (INT *)malloc(sizeof(INT)*(n+1));
    INT *newr = (INT *)malloc(sizeof(INT)*(n+1));

  
#ifdef NO_OPENMP
#pragma omp parallel for private(k)
#endif
    for(k=0;k<n;k++) {
        cq[k]= k;
    }
    group = 0;
    for(k=0;k<n;k++) {
        if ((A.IA[k+1] - A.IA[k]) > group ) group = A.IA[k+1] - A.IA[k];
    }
    *rowmax = group;

#if 0
    iavg = IA[n]/n ;	
    igold = (INT)MAX(iavg,group*0.618) +1;
    igold = group ;
#endif	

    mgl->ic = (INT *)malloc(sizeof(INT)*(group+2));
    mgl->icmap = (INT *)malloc(sizeof(INT)*(n+1));
	
    front = n-1;     
    rear = n-1;

    memset(newr, -1, sizeof(INT)*(n+1));
    memset(mgl->icmap, 0, sizeof(INT)*n);

    group=0;
    icount = 0;
    mgl->ic[0] = 0; 
    pre=0;
   
    do{
        //front = (front+1)%n;
        front ++;
        if (front == n ) front =0; // front = front < n ? front : 0 ; 
        i = cq[front];

        if(i <= pre) {
            mgl->ic[group] = icount; 
            mgl->icmap[icount] = i;
            group++; 
            icount++;
#if 0
            if ((IA[i+1]-IA[i]) > igold) 
                iend = MIN(IA[i+1], (IA[i] + igold));
	    else
#endif
            iend = IA[i+1];
            for(j= IA[i]; j< iend; j++)  newr[JA[j]] = group;
        }
        else if (newr[i] == group) { 
            //rear = (rear +1)%n;
            rear ++;  
            if (rear == n) rear = 0;
            cq[rear] = i;
        }
        else {
            mgl->icmap[icount] = i;
            icount++;
#if  0
            if ((IA[i+1] - IA[i]) > igold)  iend =MIN(IA[i+1], (IA[i] + igold));
            else
#endif 
            iend = IA[i+1];
            for(j = IA[i]; j< iend; j++)  newr[JA[j]] =  group;
        }
        pre=i;
		
//    printf("pre = %d\n",pre); 
    } while(rear != front);

//    printf("group\n"); 
    mgl->ic[group] = icount;  
    mgl->colors = group;

    free(cq);
    free(newr);
    if (theta >0 ){
    fasp_mem_free(S.IA);
    fasp_mem_free(S.JA);
    }
    *groups = group;
#endif

#if 0
    for(i=0; i < group; i++ ){
        printf( "A.color = %d A.row= %d %d\n", group, A.row, mgl->ic[i+1] - mgl->ic[i]);
    }
#endif
    return;
}


void multicolors_independent_set(AMG_data *mgl, INT gslvl)
{
#if MULTI_COLORS_ORDER
    INT i, Colors, rowmax, level, prtlvl = 3;
    REAL theta = 0.00;
    INT maxlvl = MIN(gslvl, mgl->num_levels-1);
#ifdef _OPENMP
#pragma omp parallel for private(level,rowmax,Colors) schedule(static, 1)
#endif
    for (level=0; level<maxlvl; level++) {
    //for (level=0; level<mgl->num_levels-1; level++) {

#if 0
        dCSRmat_Multicoloring(&mgl[level].A, &rowmax, &Colors);
#else
       dCSRmat_Multicoloring_Theta(&mgl[level], theta, &rowmax, &Colors);

#endif
//       if (prtlvl > 1)
//           printf("mgl[%3d].A.row = %12d rowmax = %5d rowavg = %7.2lf colors = %5d Theta = %le\n",
//           level, mgl[level].A.row, rowmax, (double)mgl[level].A.nnz/mgl[level].A.row,
//           mgl[level].colors, theta);
    }
#endif
}
#endif

