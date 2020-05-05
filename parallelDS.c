#include <stdio.h>
#include <stdlib.h>
#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"
#include "inttypes.h"
#include <mpi.h>
#include <string.h>
#include "set.h" //https://github.com/barrust/set/tree/master/src


//parameters to tune when experimenting with different blocks
#define block_size 25

#define block_h 25
#define block_w 5
// 50, 600 works well

#define FSsearch_param 3
#define DSsearch_param 3

#define pointsUnitFactor 1

#define numOfSteps 4

int min(int x, int y)
{
    return (x < y ? x : y);
}

int max(int x, int y)
{
    return (x > y ? x : y);
}

int ABS(int x,  int y)
{
    return max(x - y, y - x);
}

uint8_t *getimg(unsigned char *img, unsigned char *gray_img, size_t img_size, size_t gray_img_size, int width, int height)
{
    //function that returns a 2D array of bytes. with the byte at the (i,j) index represents the value of the pixel at that index
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
        if (((col == wb || col == wb + block_sz) && row >= hb && row <= hb + block_sz) || ((row == hb || row == hb + block_sz) && col >= wb && col <= wb + block_sz)){
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
    // int col = 0, row = 0;
    size_t size = bsize * bsize;
    uint8_t *Block = malloc(size);

    for (int i = 0; i < bsize; i++)
    {
        for (int j = 0; j < bsize; j++)
        {
            if (brow + i >= height || bcol + j >= width) { Block[i * bsize + j] = (uint8_t)(0); continue;}

            Block[i * bsize + j] = pix[(brow + i) * width + (bcol + j)];
        }
    }
    return Block;
}

// uses MSE
int evaluateBlocks(uint8_t *block1, uint8_t *block2)
{
  int diff = 0;
  for (int i = 0; i < block_size; i++)
  {
    for (int j = 0; j < block_size; j++)
    {
      diff += (int)(ABS (block1[i * block_size + j], block2[i * block_size + j])) * (ABS (block1[i * block_size + j], block2[i * block_size + j]));
    }
  }
  return diff / (block_size * block_size);
}

void sort(int size, int *arr, int *indicesX, int *indicesY)
{
    for (int i = 0; i < size; i++)
    {
        for (int j = 0; j < size - 1; j++)
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

int getMaxPixel(uint8_t *block)
{
    int mx = 0;
    for (int i = 0; i < block_size; i++){
        for (int j = 0; j < block_size; j++){
            mx = max(mx, (int)(block[i * block_size + j]));
        }
    }
    return mx;
}

//calculating peak to peak signal ratio
double PSNR(uint8_t *candidateBlock, uint8_t *originalBlock){
    int MSE = evaluateBlocks(candidateBlock, originalBlock);
    int ptp = getMaxPixel(originalBlock);
    double psnr = 20.0 * log10(ptp) - 10.0 * log10(MSE);
    return psnr;
}

int *bruteforce(uint8_t *searchGrid, uint8_t *referenceBlock, int width, int height)
{
    size_t block_si = block_size * block_size;
    uint8_t *candidateBlock = malloc(block_si);
    int best = INT_MAX;
    int *ans = (int *) malloc(3 * sizeof(int));
    printf("now exexcuting brute force...\n");
    for (int i = 0; i < FSsearch_param * block_size; i++)
    {
        for (int j = 0; j < FSsearch_param * block_size; j++)
        {
            candidateBlock = getblock(searchGrid, FSsearch_param * block_size, FSsearch_param * block_size, i, j, block_size);
            int curDiff = evaluateBlocks (candidateBlock, referenceBlock);
            // printf("brute force is now on X: %d, Y: %d, score: %d \n", i, j, curDiff);
            if (curDiff < best)
            {
              best = curDiff;
              ans[0] = i;
              ans[1] = j;
              ans[2] = curDiff;
            }
        }
    }
    printf("brute force function is done. best X: %d, Y: %d, score: %d\n", ans[0], ans[1], ans[2]);
    return ans;
}

const char *coordinatesToString (int num1, int num2){
  char *snum1 = malloc(50);
  sprintf(snum1, "%d", num1);

  char *snum2 = malloc(50);
  sprintf(snum2, "%d", num2);

  strcat (snum1, "#");
  strcat (snum1, snum2);
	printf("%s\n", snum1 );
  return snum1;
}

// int finalEXB;
// int *finalEXBcoordinates;

// here delete this function
/*
void addEXBcoordinates(int *arr, int size){
  int numOfElements = size * 2;
  int finalEXBelements = finalEXB * 2;
  for (int i = 0; i < numOfElements; i += 2){
    int found = 0;
    for (int j = 0; j < finalEXBelements; j++){
      if (arr[i] == finalEXBcoordinates[j] && arr[i + 1] == finalEXBcoordinates[j + 1]){
        found = 1;
        break;
      }
    }
    if (found == 0) {
      finalEXBcoordinates[finalEXB] = arr[i];
      finalEXBcoordinates[finalEXB + 1] = arr[i + 1];
      finalEXB += 2;
    }
  }
}
*/


int stepsCount = 0;
int localTotalEXB = 0;
int *EXBcoordinates;
int EXPsIndex = 0;

int *DSP (uint8_t *searchGrid, uint8_t *referenceBlock, int width, int height, int numOfPoints, int **points, int relativeCurCenterBlockX, int relativeCurCenterBlockY, int phase, int rank){
  size_t block_si = block_size * block_size;
  uint8_t *candidateBlock = malloc(block_si);
  int best = INT_MAX;
  int *ans;
  if (phase == 0){
    printf("proc %d.... phase 0.... current center %d, %d\n", rank, relativeCurCenterBlockX, relativeCurCenterBlockY);
    ans = (int *) malloc(9 * 3 * sizeof(int)); // array to store 9 x-y coordinates of blocks and the their 9 scores
    printf("proc %d... phase 0. created array for coordinates and scores (size 9*3)\n", rank);
    printf("proc %d... phase 0.... method addding points to initial list\n", rank);
  }
  else if (phase == 1){
    stepsCount += 1;
    printf("proc %d... phase 1... ", rank);
    if (numOfPoints == 9){
      printf("LDSP");
    }
    else{
      printf("SDSP");
    }
    printf(" %d current center %d, %d\n", stepsCount, relativeCurCenterBlockX, relativeCurCenterBlockY);
    ans = (int *) malloc(3 * sizeof(int));
  }
  for (int i = 0; i < numOfPoints; i++){
    int curPointX = relativeCurCenterBlockX + points[i][0] * pointsUnitFactor;
    int curPointY = relativeCurCenterBlockY + points[i][1] * pointsUnitFactor;

    candidateBlock = getblock(searchGrid, DSsearch_param * block_size, DSsearch_param * block_size, curPointX, curPointY, block_size);
    int curDiff = evaluateBlocks (candidateBlock, referenceBlock);

    localTotalEXB += 1; //for explored block
    EXBcoordinates[EXPsIndex] =  curPointX;
    EXBcoordinates[EXPsIndex + 1] =  curPointY;
    EXPsIndex += 2;

    if (phase == 0){
      // add coordinates and def to array
      ans[i * 3] = curPointX;
      ans[i * 3 + 1] = curPointY;
      ans[i * 3 + 2] = curDiff;
      printf("proc %d... possible point #%d. x: %d, y: %d, score: %d\n", rank, i, ans[i * 3], ans[i * 3 + 1], ans[i * 3 + 2]);
    }
    else if (phase == 1){
      printf("proc %d... phase 1... ", rank);
      if (numOfPoints == 9){
        printf("LDSP");
      }
      else{
        printf("SDSP");
      }
      printf(" %d candidate block (relative) %d, %d, %d\n", stepsCount, curPointX, curPointY, curDiff);
      if (curDiff < best)
      {
        best = curDiff;
        ans[0] = curPointX;
        ans[1] = curPointY;
        ans[2] = curDiff;
      }
    }
  }

  if (phase == 0){
    printf("proc %d... exiting phase 0...\n", rank);
    return ans;
  }

  printf("proc %d... phase 1... ", rank);
  if (numOfPoints == 9){
      printf("LDSP");
  }
  else{
    printf("SDSP");
  }
  printf(" %d center moves to %d, %d. with score %d\n\n\n\n", stepsCount, ans[0], ans[1], ans[2]);

  if ((stepsCount <= (numOfSteps - 1)) && (numOfPoints == 9 && (ans[0] != relativeCurCenterBlockX || ans[1] != relativeCurCenterBlockY))){
    ans = DSP (searchGrid, referenceBlock, width, height, numOfPoints, points, ans[0], ans[1], 1, rank);
  }

  return ans;
}

// calculate number of inspections and average among different cases
int *diamondSearch(uint8_t *searchGrid, uint8_t *referenceBlock, int width, int height,int relativeCurCenterBlockX, int relativeCurCenterBlockY, int phase, int rank)
{
  int LDSPnumOfPoints = 9;

  int **LDSPpoints;
  LDSPpoints = malloc(LDSPnumOfPoints * sizeof(int *));
  for (int i = 0; i < LDSPnumOfPoints; i++) {
    LDSPpoints[i] = malloc(2 * sizeof(int));
  }
  LDSPpoints[0][0] = 0; LDSPpoints[0][1] = 0;
  LDSPpoints[1][0] = 2; LDSPpoints[1][1] = 0;
  LDSPpoints[2][0] = -2; LDSPpoints[2][1] = 0;
  LDSPpoints[3][0] = 0; LDSPpoints[3][1] = 2;
  LDSPpoints[4][0] = 0; LDSPpoints[4][1] = -2;
  LDSPpoints[5][0] = 1; LDSPpoints[5][1] = 1;
  LDSPpoints[6][0] = -1; LDSPpoints[6][1] = 1;
  LDSPpoints[7][0] = 1; LDSPpoints[7][1] = -1;
  LDSPpoints[8][0] = -1; LDSPpoints[8][1] = -1;


  int SDSPnumOfPoints = 5;
  int **SDSPpoints;
  SDSPpoints = malloc(SDSPnumOfPoints * sizeof(int *));
  for (int i = 0; i < SDSPnumOfPoints; i++) {
    SDSPpoints[i] = malloc(2 * sizeof(int));
  }
  SDSPpoints[0][0] = 0; SDSPpoints[0][1] = 0;
  SDSPpoints[1][0] = 1; SDSPpoints[1][1] = 0;
  SDSPpoints[2][0] = -1; SDSPpoints[2][1] = 0;
  SDSPpoints[3][0] = 0; SDSPpoints[3][1] = 1;
  SDSPpoints[4][0] = 0; SDSPpoints[4][1] = -1;

  int DSsearchGridPivotCol = max(block_w - block_size * ((DSsearch_param - 1) / 2), 0);
  int DSsearchGridPivotRow = max(block_h - block_size * ((DSsearch_param - 1) / 2), 0);

  // int relativeCurCenterBlockX = ((DSsearch_param - 1) / 2) * block_size; //relative to search grid
  // int relativeCurCenterBlockY = ((DSsearch_param - 1) / 2) * block_size; //relative to search grid

  int *LDSPans;
  if (phase == 0){
    printf("proc %d... entering phase 0...\n\n", rank);

    LDSPans = (int *) malloc(18 * sizeof(int));
    LDSPans = DSP (searchGrid, referenceBlock, width, height, LDSPnumOfPoints, LDSPpoints, relativeCurCenterBlockX, relativeCurCenterBlockY, 0, rank);

    printf("\nproc %d... done phase 0... im returning 18 values\n\n", rank);

    return LDSPans;
  }
  else if (phase == 1){
    printf("proc %d... entering phase 1...\n\n", rank);

    LDSPans = (int *) malloc(3 * sizeof(int));
    printf("process %d... im inspecting this point X: %d, Y: %d\n", rank, relativeCurCenterBlockX, relativeCurCenterBlockY);
    LDSPans = DSP (searchGrid, referenceBlock, width, height, LDSPnumOfPoints, LDSPpoints, relativeCurCenterBlockX, relativeCurCenterBlockY, 1, rank);

    int *SDSPans = (int *) malloc(3 * sizeof(int));
    SDSPans = DSP (searchGrid, referenceBlock, width, height, SDSPnumOfPoints, SDSPpoints, LDSPans[0], LDSPans[1], 1, rank);

    printf("proc %d... Number of LSDP iterations = %d\n\n", rank, stepsCount - 1);

    printf("proc %d... Refernce Block location %d, %d\n", rank, DSsearchGridPivotCol + relativeCurCenterBlockX, DSsearchGridPivotRow + relativeCurCenterBlockY);
    printf("proc %d... New Block location %d, %d\n", rank, DSsearchGridPivotCol + SDSPans[0], DSsearchGridPivotRow + SDSPans[1]);

    return SDSPans;
  }
  return NULL;
}

int main(int argc, char** argv)
{
  /*
  // to do:
  // check check results of different blocks and images
  // check on multiple processors // currently works for numOfProc = 9 or less
  // comment code
  // handle memory and leaks
  // later: distribute processors on different levels and depths to get better exploration of possibilities and to allow more processors in the program
  */

  int rank, numOfProc;
  // char version [MPI_MAX_LIBRARY_VERSION_STRING];
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &numOfProc);

  EXBcoordinates = (int *)malloc(2 * FSsearch_param * FSsearch_param * block_size * block_size * sizeof(int)); // set to maximum possible // it wont reach exb of fs, but to be safe.

  //finalEXBcoordinates = (int *)malloc(2 * FSsearch_param * FSsearch_param * block_size * block_size * sizeof(int));

  // variables used by master and slaves
  int width, height, channels;
  if (rank == 0){
    printf("the MASTER is executing...\n");
  }
  else {
    printf("slave %d is executing...\n", rank);
  }

  printf("im proc %d, the size of the world is %d\n\n\n", rank, numOfProc);
  if (rank == 0){
    unsigned char *img1 = stbi_load("frame1.jpg", &width, &height, &channels, 0);
    unsigned char *img2 = stbi_load("frame2.jpg", &width, &height, &channels, 0);

    if(img1 == NULL || img2 == NULL) {
      printf("Error in loading the image\n");
      exit(1);
    }

    printf("Loaded two image with a width of %d px, a height of %d px and %d channels\n", width, height, channels);

    size_t img_size = width * height * channels;

    int gray_channels = channels == 4 ? 2 : 1;
    size_t gray_img_size = width * height * gray_channels;

    unsigned char *gray_img1 = malloc(gray_img_size);
    unsigned char *gray_img2 = malloc(gray_img_size);

    //image size 640x360
    uint8_t *pix1 = malloc(gray_img_size);
    uint8_t *pix2 = malloc(gray_img_size);

    pix1 = getimg(img1, gray_img1, img_size, gray_img_size, width, height);
    pix2 = getimg(img2, gray_img2, img_size, gray_img_size, width, height);

    printf("converted the two images to gray and stored in arrays\n");

    size_t block_si = block_size * block_size;

    uint8_t *referenceBlock = malloc(block_si);
    printf("allocated reference block size %d x %d = %d\n", block_size, block_size, block_size * block_size);

    referenceBlock = getblock(pix1, width, height, block_w, block_h, block_size);
    printf("stored reference block (upper left = %d x %d) in array\n", block_w, block_h);

    highlight_block(img1, img_size, block_w, block_h, width, height, block_size, 2);
    printf("highlighted reference block with color 2\n");



    uint8_t *FSsearchGrid = malloc(FSsearch_param * FSsearch_param * block_si);
    int FSsearchGridPivotColPoint = max(block_w - block_size * ((FSsearch_param - 1) / 2), 0);
    int FSsearchGridPivotRowPoint = max(block_h - block_size * ((FSsearch_param - 1) / 2), 0);
    FSsearchGrid = getblock(pix2, width, height, FSsearchGridPivotColPoint, FSsearchGridPivotRowPoint, block_size * FSsearch_param);
    printf("allocated space for FS searh grid (from X: %d, Y: %d)\n", FSsearchGridPivotColPoint, FSsearchGridPivotRowPoint);

    int *FSans = (int *)malloc(3 * sizeof(int));
    FSans = bruteforce(FSsearchGrid, referenceBlock, width, height);
    printf("brute force returned best X: %d, Y: %d, score: %d\n", FSans[0], FSans[1], FSans[2]);



    uint8_t *DSsearchGrid = malloc(DSsearch_param * DSsearch_param * block_si);

    int DSsearchGridPivotCol = max(block_w - block_size * ((DSsearch_param - 1) / 2), 0);
    int DSsearchGridPivotRow = max(block_h - block_size * ((DSsearch_param - 1) / 2), 0);

    printf("allocated space for DS search grid. size = %d x %d x %d = %d \n", DSsearch_param, DSsearch_param, block_size * block_size, DSsearch_param * DSsearch_param * block_size * block_size);

    DSsearchGrid = getblock(pix2, width, height, DSsearchGridPivotCol, DSsearchGridPivotRow, block_size * DSsearch_param);
    printf("stored DS search grid (upper left = %d x %d) in array\n", DSsearchGridPivotCol, DSsearchGridPivotRow);

    highlight_block(img1, img_size, DSsearchGridPivotCol, DSsearchGridPivotRow, width, height, DSsearch_param * block_size, 3);
    highlight_block(img2, img_size, DSsearchGridPivotCol, DSsearchGridPivotRow, width, height, DSsearch_param * block_size, 3);
    printf("highlighted search grid in two images\n\n");

    int relativeCurCenterBlockX = ((DSsearch_param - 1) / 2) * block_size; //relative to search grid
    int relativeCurCenterBlockY = ((DSsearch_param - 1) / 2) * block_size; //relative to search grid

    // get scores for initial 9 points to distribute to slaves
    int *DSansPossiblePoints = (int *)malloc(9 * 3 * sizeof(int)); // this should be 9 intially and 2 ba3den
    printf("allocated space for initial points identifcation. starting now... \n");
    DSansPossiblePoints = diamondSearch(DSsearchGrid, referenceBlock, width, height, relativeCurCenterBlockX, relativeCurCenterBlockY, 0, rank);

    int *initPossiblePointsX = (int *)malloc(9 * sizeof(int));
    int *initPossiblePointsY = (int *)malloc(9 * sizeof(int));
    int *initPossiblePointsScore = (int *)malloc(9 * sizeof(int));

    printf("result values returned in initial possible points.... storing in arrays\n");
    for (int i = 0; i < 9; i++){
      initPossiblePointsX[i] = DSansPossiblePoints[i * 3];
      initPossiblePointsY[i] = DSansPossiblePoints[i * 3 + 1];
      initPossiblePointsScore[i] = DSansPossiblePoints[i * 3 + 2];

      printf("x: %d, y: %d, score: %d\n", initPossiblePointsX[i], initPossiblePointsY[i], initPossiblePointsScore[i]);
    }
    printf("\n\nmoving to sorting stage\n");

    //sort according to score and keeping X and Y in order
    sort (9, initPossiblePointsScore, initPossiblePointsX, initPossiblePointsY);
    printf("now sorted ...\n");
    for (int i = 0; i < 9; i++){
      printf("x: %d, y: %d, score: %d\n", initPossiblePointsX[i], initPossiblePointsY[i], initPossiblePointsScore[i]);
    }

    printf("\n\n");

    for (int i = 1; i < numOfProc; i++){
      printf("started sending to slave %d\n", i);
      MPI_Send(referenceBlock, block_si, MPI_UINT8_T, i, 0, MPI_COMM_WORLD);
      printf("im master. sending to %d. sent referenceBlock\n", i);
      MPI_Send(DSsearchGrid, DSsearch_param * DSsearch_param * block_si, MPI_UINT8_T, i, 0, MPI_COMM_WORLD);
      printf("im master. sending to %d. sent DSsearchGrid\n", i);

      MPI_Send(&width, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
      printf("im master. sending to %d. sent width %d\n", i, width);
      MPI_Send(&height, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
      printf("im master. sending to %d. sent height %d\n", i, height);

      //send coordinates starting to proc i from sorted array
      int possibleCenterX = initPossiblePointsX[i];
      int possibleCenterY = initPossiblePointsY[i];
      MPI_Send(&possibleCenterX, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
      printf("im master. sending to %d. sent center X %d \n", i, possibleCenterX);
      MPI_Send(&possibleCenterY, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
      printf("im master. sending to %d. sent center Y %d \n", i, possibleCenterY);

      printf("ended sending to slave %d\n\n", i);
    }

    // master trying one path itself update one case yourself [0]
    int possibleCenterX = initPossiblePointsX[0];
    int possibleCenterY = initPossiblePointsY[0];
    int *DSansMaster = (int *)malloc(3 * sizeof(int));
    printf("proc %d... allocated memory for local results\n", rank);

    printf("proc %d... starting now (from %d, %d)... \n", rank, possibleCenterX, possibleCenterY);
    DSansMaster = diamondSearch(DSsearchGrid, referenceBlock, width, height, possibleCenterX, possibleCenterY, 1, rank);
    printf("proc %d... local resulted returned X:%d, Y:%d, score:%d \n", rank, DSansMaster[0], DSansMaster[1], DSansMaster[2]);

    // selecting best results results
    int *DSans = (int *)malloc(3 * sizeof(int)); // to store results from master and slaves

    DSans[0] = DSansMaster[0]; // intially set to X of master (later will be compared with results of slaves)
    DSans[1] = DSansMaster[1]; // intially set to Y of master (later will be compared with results of slaves)
    DSans[2] = DSansMaster[2]; // intially set to score of master (later will be compared with results of slaves)

    printf("proc %d... used master result to set intial DS results (will be updated as the master starts receiving results from slaves and comparing them to current vals) \n", rank);

    // exb for brute force and seq ds
    int bruteForceEXB = FSsearch_param * FSsearch_param * block_size * block_size;
    int sequentialDSEXB = localTotalEXB; // this also means that it's for the master only

    bruteForceEXB = bruteForceEXB * 1;
    sequentialDSEXB = sequentialDSEXB * 1;


    // int maxSlavesEXB = localTotalEXB; // intially set to master's exb;

    SimpleSet set;
    set_init(&set);
    printf("set init\n");


    int EXBcoordinatesElements = localTotalEXB * 2;
    for (int i = 0; i < EXBcoordinatesElements; i += 2){ //adding the elements of the master
      int coordianteX = EXBcoordinates[i];
      int coordianteY = EXBcoordinates[i + 1];

      set_add(&set, coordinatesToString(coordianteX, coordianteY));
    }

    // to receiving results from slaves
    for (int i = 1; i < numOfProc; i++){
      int *DSansSlave = (int *)malloc(3 * sizeof(int));
      MPI_Recv(DSansSlave, 3, MPI_INT, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

      printf("master... now received results from %d. results: X: %d, Y: %d, score: %d\ncurrent best values X: %d, Y: %d, score: %d\n", i, DSansSlave[0], DSansSlave[1], DSansSlave[2], DSans[0], DSans[1], DSans[2]);

      // puts result with best score in DSans
      if (DSansSlave[2] < DSans[2]){
        DSans[0] = DSansSlave[0];
        DSans[1] = DSansSlave[1];
        DSans[2] = DSansSlave[2];
        printf("proc %d... is now the current best result. X: %d, Y: %d, score: %d\n", rank, DSans[0], DSans[1], DSans[2]);
      }

      int curSlaveEXB = 0;
      MPI_Recv(&curSlaveEXB, 1, MPI_INT, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      printf("im master... receiving from slave %d... its EXB is %d.\n\n", i, curSlaveEXB);

      int *curSlaveEXBcoordinates = (int *)malloc(2 * curSlaveEXB * sizeof(int));

      MPI_Recv(curSlaveEXBcoordinates, curSlaveEXB * 2, MPI_INT, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      printf("im master... receiving from slave %d... I received array of exb\n\n", i);


      int curSlaveEXBcoordinatesElements = curSlaveEXB * 2;

      for (int j = 0; j < curSlaveEXBcoordinatesElements; j += 2){ //adding the elements of the master
        int coordianteX = curSlaveEXBcoordinates[j];
        int coordianteY = curSlaveEXBcoordinates[j + 1];

        set_add(&set, coordinatesToString(coordianteX, coordianteY));
      }
    }

    int parallelDSEXB = (int) set_length(&set);
    set_destroy(&set);

    int resultBlockX = min(DSsearchGridPivotCol + DSans[0], width - 1);
    int resultBlockY = min(DSsearchGridPivotRow + DSans[1], height - 1);
    highlight_block(img2, img_size, resultBlockX, resultBlockY, width, height, block_size, 2);
    printf("highlighted result block. upper left = %d x %d\n", resultBlockX, resultBlockY);
    stbi_write_jpg("frame_highlighted1.jpg", width, height, channels, img1, 100);
    printf("printing first frame\n");
    stbi_write_jpg("frame_highlighted2.jpg", width, height, channels, img2, 100);
    printf("printing second frame\n");




    //EXB: number of explored blocks
    printf("\n\nfinal results: the EXB of the solutions are:\nbruteForce: %d\nsequential diamond search: %d\nparallel diamond search: %d\n\n\n", bruteForceEXB, sequentialDSEXB, parallelDSEXB);



    //PSNR: measure of the accuracy/quality
    uint8_t *sequentialDSAnswer = malloc(block_si);
    sequentialDSAnswer = getblock(DSsearchGrid, DSsearch_param * block_size, DSsearch_param * block_size, DSansMaster[0], DSansMaster[1], block_size);
    double sequentialDSPSNR = PSNR(sequentialDSAnswer, referenceBlock);

    uint8_t *parallelDSAnswer = malloc(block_si);
    parallelDSAnswer = getblock(DSsearchGrid, DSsearch_param * block_size, DSsearch_param * block_size, DSans[0], DSans[1], block_size);
    double parallelDSPSNR = PSNR(parallelDSAnswer, referenceBlock);

    uint8_t *bruteForceAnswer = malloc(block_si);
    bruteForceAnswer = getblock(DSsearchGrid, DSsearch_param * block_size, DSsearch_param * block_size, FSans[0], FSans[1], block_size);
    double bruteForcePSNR = PSNR(bruteForceAnswer, referenceBlock);

    printf("\n\nfinal results: the PSNRs of the solutions are:\nbruteForce: %f\nsequential diamond search: %f\nparallel diamond search: %f\n\n", bruteForcePSNR, sequentialDSPSNR, parallelDSPSNR);











    printf("MASTER OUT! *drops mic*\n");
  }
  //slaves
  else if (rank != 0) {
    // do this for all the different points
    size_t block_si = block_size * block_size;

    uint8_t *referenceBlock = malloc(block_si);
    MPI_Recv(referenceBlock, block_si, MPI_UINT8_T, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    printf("im slave %d, i received the referenceBlock\n",rank);

    uint8_t *DSsearchGrid = malloc(DSsearch_param * DSsearch_param * block_si);
    MPI_Recv(DSsearchGrid, DSsearch_param * DSsearch_param * block_si, MPI_UINT8_T, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    printf("im slave %d, i received the DSsearchGrid\n",rank);

    MPI_Recv (&width, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    printf("im slave %d, i received the width %d\n",rank, width);

    MPI_Recv (&height, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    printf("im slave %d, i received the height %d\n",rank, height);

    int possibleCenterX;
    int possibleCenterY;

    MPI_Recv (&possibleCenterX, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    printf("im slave %d, i received the my center X %d\n",rank, possibleCenterX);

    MPI_Recv (&possibleCenterY, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    printf("im slave %d, i received the my center Y %d\n",rank, possibleCenterY);


    int *DSans = (int *)malloc(3 * sizeof(int));
    printf("proc %d... allocated memory for local results\n", rank);

    printf("proc %d... starting now (from %d, %d)... \n", rank, possibleCenterX, possibleCenterY);

    DSans = diamondSearch(DSsearchGrid, referenceBlock, width, height, possibleCenterX, possibleCenterY, 1, rank);

    printf("proc %d... local resulted returned X:%d, Y:%d, score:%d \n", rank, DSans[0], DSans[1], DSans[2]);

    // sending results to master
    printf("proc %d... now sending local results to master\n", rank);
    MPI_Send(DSans, 3, MPI_INT, 0, 0, MPI_COMM_WORLD);


    MPI_Send(&localTotalEXB, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
    printf("proc %d... sent local total EXB = %d to master \n", rank, localTotalEXB);

    MPI_Send(EXBcoordinates, localTotalEXB * 2, MPI_INT, 0, 0, MPI_COMM_WORLD);
    printf("proc %d... sent local EXB array to master \n", rank);

    printf("slave %d done. Im free at LAST!\n",rank);
  }
  MPI_Finalize();
  return 0;
}
// laklala
