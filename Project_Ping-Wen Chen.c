//
//  main.cpp
//  A
//
//  Created by Kevin Chen on 2020/2/5.
//  Copyright © 2020 Kevin Chen. All rights reserved.
//  Project 1 final

#include<stdio.h>
#include<stdlib.h>
#include<string.h>
//#include"scoringmetric.c"

int maxnum(int , int , int );
int maxcheck(int , int , int );
int maxnumLocal(int , int  ,int );
void reverse_array(int arr[], int length);

struct DP_cell{
    int score;
    char row;
    char column;
    int check;
    int lefttop;
    int left;
    int top;
    struct node *next;
    struct node *pre;
    //....
};

//struct DP_cell_re{
//    int score;
//    char row;
//    char column;
//    int check;
//    int lefttop;
//    int left;
//    int top;
//    int similar;
//};

int main(){
    FILE *input, *parameter;
    parameter = fopen("/Users/kevinchen/Desktop/WSU_Spring2020/Computational Genomics/parameter.config","r");

    int p = 0, pp = 0,sta = 1, match = 0, mismatch = 0,g = 0,h = 0;
    char par[50];

    if(parameter != NULL){
        while(fscanf(parameter,"%c", &par[p]) != EOF){
            if (par[p] == '\0' && par[p-1] == '\n') {
                break;
            }
            p++;
        }
    }
    while (pp < sizeof(par)) {
        if (par[pp] == '\n') {
            sta++;
        }
        if(par[pp] != ' ' && par[pp-1] == ' '){ // positive number
            if(sta == 1){match = atoi(&par[pp]);}
            else if(sta == 2){mismatch = atoi(&par[pp]);}
            else if(sta == 3){h = atoi(&par[pp]);}
            else if(sta == 4){g = atoi(&par[pp]);}
        }
        else if(par[pp] != ' ' && par[pp-1] != '-'){
            if(sta == 1){match = -1* atoi(&par[pp]);}
            else if(sta == 2){mismatch = -1 * atoi(&par[pp]);}
            else if(sta == 3){h = -1 * atoi(&par[pp]);}
            else if(sta == 4){g = -1 * atoi(&par[pp]);}
        }
        pp++;
    }
    printf("%d %d %d %d\n",match,mismatch,h,g);
    input = fopen("/Users/kevinchen/Desktop/WSU_Spring2020/Computational Genomics/Input1.fasta","r");
    
//    int match = 1, mismatch = -2, h = -5, g = -2; //scoring metric
    int lefttop = 0, left = 0, top = 0, flag = 0;
    int status = 1, status_title = 0, status_final = 0; //: status_t 0=title ,next換行計數器
    int i = 0, seq_count = 0, s1_c = 0, s2_c = 0,s1_ct = 0, s2_ct = 0, j = 0;
    char seq_whole[300];
    char seq1[150], seq1_T[20];
    char seq2[150], seq2_T[20];
    //讀file到統合陣列
    memset(par, ' ', sizeof(par));
    memset(seq_whole, ' ', sizeof(seq_whole));
    if(input != NULL){
        while(fscanf(input,"%c", &seq_whole[i]) != EOF){
            i++;
        }
    }
    memset(seq1, ' ', sizeof(seq1));
    memset(seq2, ' ', sizeof(seq2));
    while (seq_count<sizeof(seq_whole)) {
        //判斷是title與否
        if (seq_whole[seq_count] == '>') {status_title = 1;}
        else if(seq_whole[seq_count] == '\n'){status_title = -1;}
        //判斷第幾個seq
        if (seq_whole[seq_count] == '\n' && seq_whole[seq_count-1] == '\n') {status = 2;}
        status_final = status_title * status;
        switch (status_final) {
            case 1: //sequence 1 information
                if (seq_whole[seq_count] != '\n')
                {strncpy(&seq1_T[s1_ct], &seq_whole[seq_count],1);s1_ct++;}
                break;
            case -1: // sequence 1
                if (seq_whole[seq_count] != '\n')
                {
                    if (seq_whole[seq_count] == 'A' || seq_whole[seq_count] == 'C' || seq_whole[seq_count] == 'G' || seq_whole[seq_count] == 'T'||seq_whole[seq_count] == 'a' || seq_whole[seq_count] == 'c' || seq_whole[seq_count] == 'g' || seq_whole[seq_count] == 't') {
                        strncpy(&seq1[s1_c] ,&seq_whole[seq_count],1);
                        s1_c++;
                    }
                }
                break;
            case 2: //sequence 2 information
                if (seq_whole[seq_count] != '\n')
                {strncpy(&seq2_T[s2_ct], &seq_whole[seq_count],1);s2_ct++;}
                break;
            case -2: //sequence 2
                if (seq_whole[seq_count] != '\n')
                {
                    if (seq_whole[seq_count] == 'A' || seq_whole[seq_count] == 'C' || seq_whole[seq_count] == 'G' || seq_whole[seq_count] == 'T'||seq_whole[seq_count] == 'a' || seq_whole[seq_count] == 'c' || seq_whole[seq_count] == 'g' || seq_whole[seq_count] == 't') {
                        strncpy(&seq2[s2_c] ,&seq_whole[seq_count],1);
                        s2_c++;
                    }
                }
                break;
            default:
                break;
        }
        seq_count++;
    }
    printf("%s\n",seq_whole);
    printf("%s\n",seq1);
    printf("%s\n",seq2);

    
    // struct operation----------------------------------------------------
    struct DP_cell cell[s2_c+1][s1_c+1];
    
    cell[0][0].score = 0;
    for (i=1; i<s2_c+1; i++) {
        cell[i][0].score = i*g + h;
    }
    for (j=1; j<s1_c+1; j++) {
        cell[0][j].score = j*g + h;
    }
    for (i=1; i<s2_c+1; i++) {
        for (j=1; j<s1_c+1; j++) {
            if(i!=0 && j!=0){
                strncpy(&cell[i][j].row, &seq2[i-1], 1);
                strncpy(&cell[i][j].column, &seq1[j-1], 1);
            }
        }
    }

    // calculate metric table----------------------------------------------
    // calculate metric table----------------------------------------------
    int choice = 0;
    printf("1.Gloabl alignment\n2.Local alignment\n");
    scanf("%d",&choice);
    switch (choice) {
        case 1:
            for(i=1;i<s2_c+1;i++){
                for(j=1;j<s1_c+1;j++){
                    if(cell[i][j].row == cell[i][j].column){ //match
                        lefttop = cell[i-1][j-1].score + match; //lefttop + match
                        //left score check
                        if(cell[i][j-1].check == 1){ left = cell[i][j-1].score + g + h;}
                        else if(cell[i][j-1].check == 2){ left = cell[i][j-1].score + g;}
                        else{ left = cell[i][j-1].score + h + g;}
                        //top score check
                        if(cell[i-1][j].check == 1){ top = cell[i-1][j].score + g + h;}
                        else if(cell[i-1][j].check == 2){ top = cell[i-1][j].score + g + h;}
                        else{ top = cell[i-1][j].score + g;}
                        // find max of three value, assign 1 row 1 column score by adding right from top left.
                        cell[i][j].lefttop = lefttop;
                        cell[i][j].left = left;
                        cell[i][j].top = top;
                        if(i==1 || j==1){ //Global condition
                            cell[i][j].score = cell[i-1][j-1].score + match;
                            cell[i][j].check = 1;
                        }
                        else{//find max value
                            cell[i][j].score = maxnum(lefttop, left ,top);
                            cell[i][j].check = maxcheck(lefttop, left ,top); // 1 lefttop, 2 left, 3 top
                        }
                    }
                    else{//mismatch
                        lefttop = cell[i-1][j-1].score + mismatch;// lefttop + mismatch
                        //left score check
                        if(cell[i][j-1].check == 1){ left = cell[i][j-1].score + g + h;}
                        else if(cell[i][j-1].check == 2){ left = cell[i][j-1].score + g;}
                        else{ left = cell[i][j-1].score + h + g;}
                        //top score check
                        if(cell[i-1][j].check == 1){ top = cell[i-1][j].score + g + h;}
                        else if(cell[i-1][j].check == 2){ top = cell[i-1][j].score + g + h;}
                        else{ top = cell[i-1][j].score + g;}
                        // find max of three value, assign 1 row 1 column score by adding right from top left.
                        cell[i][j].lefttop = lefttop;
                        cell[i][j].left = left;
                        cell[i][j].top = top;
                        if(i==1 || j==1){ // Global condition
                            cell[i][j].score = cell[i-1][j-1].score + match;
                            cell[i][j].check = 1;
                        }
                        else{ //find max value
                            cell[i][j].score = maxnum(lefttop, left ,top);
                            cell[i][j].check = maxcheck(lefttop, left ,top); // 1 lefttop, 2 left, 3 top
                        }
                    }
                }
            }
            break;
        case 2:
            for(i=1;i<s2_c+1;i++){
                for(j=1;j<s1_c+1;j++){
                    if(cell[i][j].row == cell[i][j].column){ //match
                        lefttop = cell[i-1][j-1].score + match;
                        if(cell[i][j-1].check == 1){ left = cell[i][j-1].score + g + h;}
                        else if(cell[i][j-1].check == 2){ left = cell[i][j-1].score + g;}
                        else{ left = cell[i][j-1].score + h + g;}
                        if(cell[i-1][j].check == 1){ top = cell[i-1][j].score + g + h;}
                        else if(cell[i-1][j].check == 2){ top = cell[i-1][j].score + g + h;}
                        else{ top = cell[i-1][j].score + g;}
                        // find max of three value, assign 1 row 1 column score by adding right from top left.
                        cell[i][j].lefttop = lefttop;
                        cell[i][j].left = left;
                        cell[i][j].top = top;
                        cell[i][j].score = maxnumLocal(lefttop, left ,top);
                        cell[i][j].check = maxcheck(lefttop, left ,top); // 1 lefttop, 2 left, 3 top
                    }
                    else{//mismatch
                    lefttop = cell[i-1][j-1].score + mismatch;
                        if(cell[i][j-1].check == 1){ left = cell[i][j-1].score + g + h;}
                        else if(cell[i][j-1].check == 2){ left = cell[i][j-1].score + g;}
                        else{ left = cell[i][j-1].score + h + g;}
                        if(cell[i-1][j].check == 1){ top = cell[i-1][j].score + g + h;}
                        else if(cell[i-1][j].check == 2){ top = cell[i-1][j].score + g + h;}
                        else{ top = cell[i-1][j].score + g;}
                        // find max of three value, assign 1 row 1 column score by adding right from top left.
                        cell[i][j].lefttop = lefttop;
                        cell[i][j].left = left;
                        cell[i][j].top = top;
                        cell[i][j].score = maxnumLocal(lefttop, left ,top);
                        cell[i][j].check = maxcheck(lefttop, left ,top); // 1 lefttop, 2 left, 3 top
                    }
                }
            }
            break;
        default:
            break;
    }

    //backtracking------------------------------------------------------
    //backtracking------------------------------------------------------
    
    struct DP_cell st[500];
    int count = 0, sc = 0;
    i = s2_c; j = s1_c;
    switch (choice) {
        case 1:
            while (flag == 0) {
                if (count == 0) {
                    st[count].score = cell[i][j].score;
                    strncpy(&st[count].row, &cell[i][j].row, 1);
                    strncpy(&st[count].column, &cell[i][j].column, 1);
                    st[count].check = cell[i][j].check;
                    st[count].lefttop = cell[i][j].lefttop;
                    st[count].left = cell[i][j].left;
                    st[count].top = cell[i][j].top;
                    count++;
                }
                else{
                    if (st[count-1].check == 1) {
                        i--;
                        j--;
                        st[count].score = cell[i][j].score;
                        strncpy(&st[count].row, &cell[i][j].row, 1);
                        strncpy(&st[count].column, &cell[i][j].column, 1);
                        st[count].check = cell[i][j].check;
                        st[count].lefttop = cell[i][j].lefttop;
                        st[count].left = cell[i][j].left;
                        st[count].top = cell[i][j].top;
                        count++;
                    }
                    else if (st[count-1].check == 2) {
                        j--;
                        st[count].score = cell[i][j].score;
                        strncpy(&st[count].row, &cell[i][j].row, 1);
                        strncpy(&st[count].column, &cell[i][j].column, 1);
                        st[count].check = cell[i][j].check;
                        st[count].lefttop = cell[i][j].lefttop;
                        st[count].left = cell[i][j].left;
                        st[count].top = cell[i][j].top;
                        count++;
                    }
                    else if (st[count-1].check == 3) {
                        i--;
                        st[count].score = cell[i][j].score;
                        strncpy(&st[count].row, &cell[i][j].row, 1);
                        strncpy(&st[count].column, &cell[i][j].column, 1);
                        st[count].check = cell[i][j].check;
                        st[count].lefttop = cell[i][j].lefttop;
                        st[count].left = cell[i][j].left;
                        st[count].top = cell[i][j].top;
                        count++;
                    }
                    if (i == 0 && j == 0) {
                        flag = 1;
                    }
                }
            }
            break;
        case 2:
            while (flag == 0) {
                if (count == 0) {
                    st[count].score = cell[i][j].score;
                    strncpy(&st[count].row, &cell[i][j].row, 1);
                    strncpy(&st[count].column, &cell[i][j].column, 1);
                    st[count].check = cell[i][j].check;
                    st[count].lefttop = cell[i][j].lefttop;
                    st[count].left = cell[i][j].left;
                    st[count].top = cell[i][j].top;
                    count++;
                }
                else{
                    if (st[count-1].check == 1) {
                        i--;
                        j--;
                        st[count].score = cell[i][j].score;
                        strncpy(&st[count].row, &cell[i][j].row, 1);
                        strncpy(&st[count].column, &cell[i][j].column, 1);
                        st[count].check = cell[i][j].check;
                        st[count].lefttop = cell[i][j].lefttop;
                        st[count].left = cell[i][j].left;
                        st[count].top = cell[i][j].top;
                        count++;
                    }
                    else if (st[count-1].check == 2) {
                        j--;
                        st[count].score = cell[i][j].score;
                        strncpy(&st[count].row, &cell[i][j].row, 1);
                        strncpy(&st[count].column, &cell[i][j].column, 1);
                        st[count].check = cell[i][j].check;
                        st[count].lefttop = cell[i][j].lefttop;
                        st[count].left = cell[i][j].left;
                        st[count].top = cell[i][j].top;
                        count++;
                    }
                    else if (st[count-1].check == 3) {
                        i--;
                        st[count].score = cell[i][j].score;
                        strncpy(&st[count].row, &cell[i][j].row, 1);
                        strncpy(&st[count].column, &cell[i][j].column, 1);
                        st[count].check = cell[i][j].check;
                        st[count].lefttop = cell[i][j].lefttop;
                        st[count].left = cell[i][j].left;
                        st[count].top = cell[i][j].top;
                        count++;
                    }
                }
                sc = cell[i][j].score;
                if (sc == 0) {
                    flag = 1;
                }
            }
            break;
        default:
            break;
    }

    // reverse st -----------------------------------------------------------------------
    // reverse st -----------------------------------------------------------------------
    
    int temp = count-1;
    struct DP_cell st1[count];
    for (i = 0; i < count+1; i++) {
        st1[i].check = st[temp].check;
        st1[i].score = st[temp].score;
        st1[i].lefttop = st[temp].lefttop;
        st1[i].left = st[temp].left;
        st1[i].top = st[temp].top;
        strncpy(&st1[i].row, &st[temp].row, 1);
        strncpy(&st1[i].column, &st[temp].column, 1);
        temp--;
    }
    for (i=0; i<count; i++) {
        printf("%d %c %c  %d\n",st1[i].score,st1[i].row,st1[i].column,st1[i].check);
    }
    int optscr = st1[count-1].score;

    // print result ---------------------------------------------------------------------
    // print result ---------------------------------------------------------------------
    
    int og_c = 0, g_c = 0, m = 0, mm = 0;
    int k = 0, pri = 0, flag1 = 0;
    int stage_count = 0, stage_count1 = 0, stage_count2 = 0;
    int stage = 1;
    int s2_count = 0, s1_count = 0;
    float i_percent, g_percent;
    char r, co;
    printf("Scores: match = %d, mismatch = %d, h = %d, g = %d\n",match,mismatch,h,g);
    printf("Sequence 1, length = %d characters\n",s1_c);
    printf("Sequence 2, length = %d characters\n\n",s2_c);
    pri = count / 60; i = 1;j = 1;k = 1;
    while (flag1 == 0) {
        if (stage == 1) {
            printf("%c",st1[i].column);
            if(st1[i].column != '\0'){s1_count++;}
            if(stage_count == 59){stage = 2; printf("\t%d\n",s1_count); stage_count = -1;}
            stage_count++; i++;
        }
        else if(stage == 2){
            strncpy(&r, &st1[j].row, 1);
            strncpy(&co, &st1[j].column, 1);
            if (st1[j].check == 2 || st1[j].row == '\0') {printf(" ");}
            else{
                if (co == r) {printf("|");m++;}
                else{printf(" ");mm++;}
            }
            if (stage_count1 == 59) { printf("\n");stage = 3; stage_count1 = -1;} // end of this row set stage to 3
            stage_count1++;
            j++;
        }
        else{
            if (st1[k].check == 2) {
                printf("-");g_c++;
                if (st1[k-1].check != 2) {
                    og_c++;
                }
            }
            else{
                if (st1[k].row != '\0') {
                    printf("%c",st1[k].row); //stage 3 print char and -
                    s2_count++;
                }
            }
            if (k == count-1) {printf("\t%d\n",s2_count); flag1 = 1;}
            else{
                if (stage_count2 == 59) {printf("\t%d\n",s2_count); printf("\n"); stage_count2 = -1;stage = 1;} // end of stage 3 printf next line
                stage_count2++;
                k++;
            }
        }//end of print break loop
    }
    printf("\nReport: \nGlobal optimal score = %d\nNumbers of matches: matches = %d, mismatches = %d, gaps = %d, opening gaps = %d\n",optscr,m,mm-1,g_c,og_c);
    if (s1_c>s2_c) {
        i_percent = ((float)m/s1_c)*100;
        g_percent = ((float)g_c/s1_c)*100;
    }
    else{
        i_percent = ((float)m/s2_c)*100;
        g_percent = ((float)g_c/s2_c)*100;
    }
    printf("Identities = %d/%d (%.f%%), Gaps = %d/%d (%.f%%)\n",m,s1_c,i_percent,g_c,s1_c,g_percent);
    
    fclose(input);
//    fclose(output);
    return 0;
}

int maxnum(int i, int j ,int k)
{
    int c;
    if (i > j && i > k) {
        c = i;
    }
    else if(j > k){
        c = j;
    }
    else{
        c = k;
    }
    return c;
}

int maxcheck(int i, int j, int k)
{
    int c;

    if (i > j && i > k) {
        c = 1;
    }
    else if(j > k){
        c = 2;
    }
    else{
        c = 3;
    }
    return c;
}
int maxnumLocal(int i, int j ,int k)
{
    int c;
    if (i > j && i > k) {
        if (i > 0) {
            c = i;
        }
        else{c = 0;}
    }
    else if(j > k){
        if (j > 0) {
            c = j;
        }
        else{c = 0;}
    }
    else{
        if (k > 0) {
            c = k;
        }
        else{c = 0;}
    }
    return c;
}

void reverse_array(int arr[], int length)
{
        int i, temp;
        for (i=0; i<length/2; i++)
        {
                 temp = arr[i];
                 arr[i] = arr[length-i-1];
                 arr[length-i-1] = temp;
        }
}

//while(seq_count < sizeof(seq_whole)){
//        if (seq_whole[seq_count] == '2') {
//            status = 1;
//            seq_count += 1;
//        }
//        if(seq_whole[seq_count] != '\n' && seq_whole[seq_count] != '2' && seq_whole[seq_count] != '>' && seq_whole[seq_count] != 's' && seq_whole[seq_count] != '1' && status == 0){
//            seq1[row1][s1_c] = seq_whole[seq_count];
//            printf("%c",seq1[row1][s1_c]);
//            s1_c++;
//            if(s1_c%60 == 0){row1++;printf("\n");}
//        }
//        else if(seq_whole[seq_count] != '2' && seq_whole[seq_count] != '>' && seq_whole[seq_count] != 's' && seq_whole[seq_count] != '1' && status == 1){
//            //("%c\n",seq_whole[seq_count]);
//            if (seq_whole[seq_count] == '\n'){
//                if (s2_c%60 != 0) {
//                    in_count2 = s2_c%60;
//                    for(i=0; i<60 - in_count2;i++){
//                        seq2[row2][s2_c] = ' ';
//                        s2_c++;
//                    }
//                    row2++;
//                }
//                else{
//                    //seq_count++;
//                    seq2[row2][s2_c] = seq_whole[seq_count];
//                }
//            }
//            else{
//                seq2[row2][s2_c] = seq_whole[seq_count];
//                printf("%c",seq2[row2][s2_c]);
//                s2_c++;
//                if(s2_c%60 == 0){row2++;printf("\n");}
//            }
//        }
//        seq_count++;
//    }
//cout << 55;
//System.out.println(55);
//print(55)
//console.log(55)
