#include <stdio.h>
#include <stdlib.h>
#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"
#include "inttypes.h"

//parameters to tune when experimenting with different blocks
#define block_size 25

#define block_h 50
#define block_w 600

// 50, 600 works well

#define FSsearch_param 3
#define DSsearch_param 3

#define pointsUnitFactor 1

#define numOfSteps 200

int min(int x , int y)
{
    return (x < y ? x : y);
}

int max(int x , int y)
{
    return (x > y ? x : y);
}

int ABS(int x,  int y)
{
    return max(x - y , y - x);
}

uint8_t *getimg(unsigned char *img, unsigned char *gray_img, size_t img_size, size_t gray_img_size, int width, int height)
{
    //function that return a 2D array of bytes where the byte at the (i,j) index represents the value of the pixel at that index
    //also the function will turn the image into gray scaled
    int col = 0,row = 0;
    size_t size = height * width;
    uint8_t *pix = malloc(size);
    for(unsigned char *p = img, *pg = gray_img; p != img + img_size; p += 3, pg += 1)
    {
        *pg = (uint8_t)((*p + *(p + 1) + *(p + 2)) / 3.0);
        pix[row * width + col] = (uint8_t) *pg;

        col++;
        if (col >= width) {col = 0; row++;}
    }
    return pix;
}


void highlight_block(unsigned char *curimg, size_t img_size, int wb, int hb, int width, int height, int block_sz, int color)
{
    int col = 0, row = 0;
    for(unsigned char *p = curimg; p != curimg + img_size; p += 3)
    {
        if ( ((col == wb || col == wb + block_sz) && row >= hb && row <= hb + block_sz ) || ((row == hb || row == hb + block_sz) && col >= wb && col <= wb + block_sz ) ){
          if (color == 1){
            *(p) = (uint8_t)((*p + *(p + 1) + *(p + 2)) /3.0);
            *(p+1) = (uint8_t)(0);
            *(p+2) = (uint8_t)(0);
          }
          if (color == 2){
            *(p) = (uint8_t)(0);
            *(p+1) = (uint8_t)((*p + *(p + 1) + *(p + 2)) /3.0);
            *(p+2) = (uint8_t)(0);
          }
          if (color == 3){
            *(p) = (uint8_t)(0);
            *(p+1) = (uint8_t)(0);
            *(p+2) = (uint8_t)((*p + *(p + 1) + *(p + 2)) /3.0);
          }
        }

        col++;
        if (col >= width) {col = 0; row++;}
    }
}

//only used with gray scaled imgs
//make sure to check boundaries when passing parameters bcol, brow
// takes upper left cornor and returns block as an array
  uint8_t *getblock(uint8_t *pix, int width, int height, int bcol, int brow, int bsize)
{
    int col = 0, row = 0;
    size_t size = bsize * bsize;
    uint8_t *Block = malloc(size);

    for (int i = 0; i < bsize; i++)
    {
        for (int j = 0; j < bsize; j++)
        {
            if (brow + i >= height || bcol + j >= width ) { Block[i * bsize + j] = (uint8_t)(0); continue;}

            Block[i * bsize + j] = pix[(brow + i) * width + (bcol + j)];
        }
    }
    return Block;
}

int evaluateBlocks(uint8_t *block1 , uint8_t *block2)
{
    int diff = 0;
    for (int i = 0; i < block_size; i++)
    {
        for (int j = 0; j < block_size; j++)
        {
            diff += (int)(ABS (block1[i * block_size + j], block2[i * block_size + j]) );
        }
    }
    return diff;
}

void sort(int size, int *arr, int *indicesX, int *indicesY)
{
    for (int i = 0; i < size; i++)
    {
        for (int j = 0; j < size - 1 ; j++)
        {
            if (arr[j] > arr[j + 1])
            {
                int temp = arr[j];
                arr[j] = arr[j + 1];
                arr[j + 1] = temp;

                temp = indicesX[j];
                indicesX[j] = indicesX[j + 1];
                indicesX[j + 1] = temp;

                temp = indicesY[j];
                indicesY[j] = indicesY[j + 1];
                indicesY[j + 1] = temp;
            }
        }
    }
}


int *bruteforce(uint8_t *searchGrid, uint8_t *referenceBlock, int width, int height)
{
    size_t block_si = block_size * block_size;
    uint8_t *candidateBlock = malloc(block_si);
    int best = INT_MAX;
    int *ans = (int *) malloc(2 * sizeof(int));
    for (int i = 0; i < FSsearch_param * block_size; i++)
    {
        for (int j = 0; j < FSsearch_param * block_size; j++)
        {
            candidateBlock = getblock(searchGrid, FSsearch_param * block_size, FSsearch_param * block_size, j, i, block_size);
            int curDiff = evaluateBlocks (candidateBlock, referenceBlock);
            if (curDiff < best)
            {
                ans[0] = i;
                ans[1] = j;
                best = curDiff;
            }
        }
    }
    // printf("%d %d\n" , ans[0] , ans[1]);
    return ans;
}


int counter = 0;

int *DSP (uint8_t *searchGrid, uint8_t *referenceBlock, int width, int height, int numOfPoints, int **points, int relativeCurCenterBlockX, int relativeCurCenterBlockY, int phase){
  size_t block_si = block_size * block_size;
  uint8_t *candidateBlock = malloc(block_si);
  int best = INT_MAX;
  int *ans;
  if (phase == 0){
    printf(" phase 0....  current center %d, %d\n",  relativeCurCenterBlockX , relativeCurCenterBlockY);
    printf("phase 0. created array for coordinates and SADs (size 9*3)\n");
    ans = (int *) malloc(9 * 3 * sizeof(int)); // array to store 9 x-y coordinates of blocks and the their 9 SADs
    printf("phase 0.... method addding point to initial list\n");
  }
  else if (phase == 1){
    counter += 1;
    printf("phase 1... ");
    if (numOfPoints == 9){
      printf("LDSP");
    }
    else{
      printf("SDSP");
    }
    printf(" %d current center %d, %d\n" , counter,  relativeCurCenterBlockX , relativeCurCenterBlockY);
    ans = (int *) malloc(2 * sizeof(int));
  }
  for (int i = 0; i < numOfPoints; i++){
    int curPointX = relativeCurCenterBlockX + points[i][0] * pointsUnitFactor;
    int curPointY = relativeCurCenterBlockY + points[i][1] * pointsUnitFactor;

    candidateBlock = getblock(searchGrid, DSsearch_param * block_size, DSsearch_param * block_size, curPointX, curPointY, block_size);
    int curDiff = evaluateBlocks (candidateBlock, referenceBlock);
    if (phase == 0){
      // add coordinates and def to array
      ans[i * 3] = curPointX;
      ans[i * 3 + 1] = curPointY;
      ans[i * 3 + 2] = curDiff;
      printf("x: %d, y: %d, score: %d\n", ans[i * 3], ans[i * 3 + 1], ans[i * 3 + 2]);
    }
    else if (phase == 1){
      printf("candidate block (relative) %d, %d, %d\n", curPointX, curPointY, curDiff);
      if (curDiff < best)
      {
        ans[0] = curPointX;
        ans[1] = curPointY;
        best = curDiff;
      }
    }

  }


  if (phase == 0){
    printf("exiting phase 0...\n");
    return ans;
  }

  printf("phase 1... ");
  if (numOfPoints == 9){
      printf("LDSP");
  }
  else{
    printf("SDSP");
  }
  printf(" %d center moves to %d, %d\n\n\n\n", counter, ans[0], ans[1]);


  if ((counter <= (numOfSteps - 1)) && (numOfPoints == 9 && (ans[0] != relativeCurCenterBlockX || ans[1] != relativeCurCenterBlockY))){
    ans = DSP (searchGrid, referenceBlock, width, height, numOfPoints, points, ans[0], ans[1], 1);
  }
  // ans[0] = relativeCurCenterBlockX;
  // ans[1] = relativeCurCenterBlockY;

  return ans;
}

// calculate number of inspections and average among different cases
int *diamondSearch(uint8_t *searchGrid, uint8_t *referenceBlock, int width, int height,int relativeCurCenterBlockX, int relativeCurCenterBlockY, int phase)
{
  int LDSPnumOfPoints = 9;

  int **LDSPpoints;
  LDSPpoints =  malloc(LDSPnumOfPoints * sizeof(int *));
  for (int i = 0; i < LDSPnumOfPoints; i++) {
    LDSPpoints[i] =  malloc(2 * sizeof(int));
  }
  LDSPpoints[0][0] =  0; LDSPpoints[0][1] = 0;
  LDSPpoints[1][0] =  2; LDSPpoints[1][1] = 0;
  LDSPpoints[2][0] =  -2; LDSPpoints[2][1] = 0;
  LDSPpoints[3][0] =  0; LDSPpoints[3][1] = 2;
  LDSPpoints[4][0] =  0; LDSPpoints[4][1] = -2;
  LDSPpoints[5][0] =  1; LDSPpoints[5][1] = 1;
  LDSPpoints[6][0] =  -1; LDSPpoints[6][1] = 1;
  LDSPpoints[7][0] =  1; LDSPpoints[7][1] = -1;
  LDSPpoints[8][0] =  -1; LDSPpoints[8][1] = -1;


  int SDSPnumOfPoints = 5;
  int **SDSPpoints;
  SDSPpoints =  malloc(SDSPnumOfPoints * sizeof(int *));
  for (int i = 0; i < SDSPnumOfPoints; i++) {
    SDSPpoints[i] =  malloc(2 * sizeof(int));
  }
  SDSPpoints[0][0] =  0; SDSPpoints[0][1] = 0;
  SDSPpoints[1][0] =  1; SDSPpoints[1][1] = 0;
  SDSPpoints[2][0] =  -1; SDSPpoints[2][1] = 0;
  SDSPpoints[3][0] =  0; SDSPpoints[3][1] = 1;
  SDSPpoints[4][0] =  0; SDSPpoints[4][1] = -1;

  int DSsearchGridPivotCol = max(block_w - block_size * ((DSsearch_param - 1) / 2), 0);
  int DSsearchGridPivotRow = max(block_h - block_size * ((DSsearch_param - 1) / 2), 0);

  // int relativeCurCenterBlockX = ((DSsearch_param - 1) / 2) * block_size;  //relative to search grid
  // int relativeCurCenterBlockY = ((DSsearch_param - 1) / 2) * block_size;  //relative to search grid

  int *LDSPans;
  if (phase == 0){
    printf("entering phase 0...\n\n");

    LDSPans = (int *) malloc(18 * sizeof(int));
    LDSPans = DSP (searchGrid, referenceBlock, width, height, LDSPnumOfPoints, LDSPpoints, relativeCurCenterBlockX, relativeCurCenterBlockY, 0);

    printf("\n\ndone phase 0... returned 18 values\n\n");

    return LDSPans;
  }
  else if (phase == 1){
    printf("entering phase 1...\n\n");

    LDSPans = (int *) malloc(2 * sizeof(int));
    LDSPans = DSP (searchGrid, referenceBlock, width, height, LDSPnumOfPoints, LDSPpoints, relativeCurCenterBlockX, relativeCurCenterBlockY, 1);

    int *SDSPans = (int *) malloc(2 * sizeof(int));
    SDSPans = DSP (searchGrid, referenceBlock, width, height, SDSPnumOfPoints, SDSPpoints, LDSPans[0], LDSPans[1], 1);

    printf("Number of LSDP iterations = %d\n\n", counter - 1);

    printf("Refernce Block location %d, %d\n", DSsearchGridPivotCol + relativeCurCenterBlockX, DSsearchGridPivotRow + relativeCurCenterBlockY);
    printf("New Block location %d, %d\n", DSsearchGridPivotCol + SDSPans[0], DSsearchGridPivotRow + SDSPans[1]);

    return SDSPans;
  }
}

int main(int argc, char** argv`)
{
  int rank, numOfProc;
  char version [MPI_MAX_LIBRARY_VERSION_STRING];
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &numOfProc);

  if (rank == 0){
    int width, height, channels;
    unsigned char *img1 = stbi_load("frame1.jpg", &width, &height, &channels, 0);
    unsigned char *img2 = stbi_load("frame2.jpg", &width, &height, &channels, 0);

    printf("Loaded image with a width of %d px, a height of %d px and %d channels\n", width, height, channels);

    if(img1 == NULL || img2 == NULL) {
      printf("Error in loading the image\n");
      exit(1);
    }

    size_t img_size = width * height * channels;
    int gray_channels = channels == 4 ? 2 : 1;

    size_t gray_img_size = width * height * gray_channels;
    unsigned char *gray_img1 = malloc(gray_img_size);
    unsigned char *gray_img2 = malloc(gray_img_size);

    //image size 640x360
    uint8_t *pix1 = malloc(gray_img_size);
    uint8_t *pix2 = malloc(gray_img_size);
    size_t block_si = block_size * block_size;
    uint8_t *referenceBlock = malloc(block_si);

    uint8_t *FSsearchGrid = malloc(FSsearch_param * FSsearch_param * block_si);
    uint8_t *DSsearchGrid = malloc(DSsearch_param * DSsearch_param * block_si);

    pix1 = getimg(img1 , gray_img1 , img_size , gray_img_size , width , height);
    pix2 = getimg(img2, gray_img2, img_size, gray_img_size, width, height);

    printf("converted the two images to gray and stored in arrays\n");

    referenceBlock = getblock(pix1, width, height, block_w, block_h, block_size);
    printf("stored reference block in indep array\n");
    highlight_block(img1 , img_size , block_w , block_h , width , height , block_size, 2);
    printf("highlighted reference block\n");


    //int FSsearchGridPivotColPoint = max(block_w - block_size * ((FSsearch_param - 1) / 2), 0);
    //int FSsearchGridPivotRowPoint = max(block_h - block_size * ((FSsearch_param - 1) / 2), 0);
    // FSsearchGrid = getblock(pix2, width, height, FSsearchGridPivotColPoint, FSsearchGridPivotRowPoint, block_size * FSsearch_param);

    // int *FSans = (int *)malloc(2 * sizeof(int));
    // FSans = bruteforce(FSsearchGrid, referenceBlock, width, height);


    int DSsearchGridPivotCol = max(block_w - block_size * ((DSsearch_param - 1) / 2), 0);
    int DSsearchGridPivotRow = max(block_h - block_size * ((DSsearch_param - 1) / 2), 0);
    DSsearchGrid = getblock(pix2, width, height, DSsearchGridPivotCol, DSsearchGridPivotRow, block_size * DSsearch_param);

    printf("stored DS search grid in indep array\n");
    highlight_block(img1 , img_size , DSsearchGridPivotCol , DSsearchGridPivotRow , width , height , DSsearch_param *  block_size, 3);
    highlight_block(img2 , img_size , DSsearchGridPivotCol , DSsearchGridPivotRow , width , height , DSsearch_param *  block_size, 3);
    printf("highlighted search grid in two images\n\n");

    int relativeCurCenterBlockX = ((DSsearch_param - 1) / 2) * block_size;  //relative to search grid
    int relativeCurCenterBlockY = ((DSsearch_param - 1) / 2) * block_size;  //relative to search grid


    // for master
    int *DSans = (int *)malloc(9 * 3 * sizeof(int)); // this should be 9 intially and 2 ba3den
    DSans = diamondSearch(DSsearchGrid, referenceBlock, width, height, relativeCurCenterBlockX, relativeCurCenterBlockY, 0);

    int *initPossiblePointsX = (int *)malloc(9 * sizeof(int));
    int *initPossiblePointsY = (int *)malloc(9 * sizeof(int));
    int *initPossiblePointsScore = (int *)malloc(9 * sizeof(int));

    printf("result values in DS ans.... storing in indep arrays\n");
    for (int i = 0; i < 9; i++){
      initPossiblePointsX[i] = DSans[i * 3];
      initPossiblePointsY[i] = DSans[i * 3 + 1];
      initPossiblePointsScore[i] = DSans[i * 3 + 2];

      printf("x: %d, y: %d, score: %d\n", initPossiblePointsX[i], initPossiblePointsY[i], initPossiblePointsScore[i]);
    }
    printf("\n\nmoving to sorting stage\n");
    //sort according to score and keeping X and Y in order
    sort (9, initPossiblePointsScore, initPossiblePointsX, initPossiblePointsY);
    printf("now sorted ...\n");
    for (int i = 0; i < 9; i++){
      printf("x: %d, y: %d, score: %d\n", initPossiblePointsX[i], initPossiblePointsY[i], initPossiblePointsScore[i]);
    }

    highlight_block(img2 , img_size , min(DSsearchGridPivotCol + DSans[0] , width-1) ,  min(DSsearchGridPivotRow + DSans[1] , height - 1) , width , height , block_size, 2);

    stbi_write_jpg("frame_highlighted1.jpg", width, height, channels, img1, 100);
    stbi_write_jpg("frame_highlighted2.jpg", width, height, channels, img2, 100);

  }
  //slaves
  else if (rank != 0) {
    int *DSans = (int *)malloc(2 * sizeof(int));
    DSans = diamondSearch(DSsearchGrid, referenceBlock, width, height, relativeCurCenterBlockX, relativeCurCenterBlockY, 1);
  }
  MPI_Finalize();
  return 0;
}
