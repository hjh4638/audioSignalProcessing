#include <stdio.h>
#define _USE_MATH_DEFINES
#include <math.h>

float x[10000] = {0,};
int COUNT = 0;
void getAudio(){
	FILE	*fi;
	float	input;
	short	data;

	fi = fopen("2.raw", "rb");
	for( ; ; ) {
		if(fread(&data, 2, 1, fi) == NULL)
			break;

		input = (float)data;  // DFT에서는 array로 저장
		x[COUNT] = input;
		COUNT++;
	}
	fclose(fi);
}
int main(){
	float Xr[10000] = {0,};
	float Xi[10000] = {0,};
	float X[10000] = {0,};
	
	int N = 1000;

	FILE* fout;
	fout = fopen("2out.txt", "w");

	getAudio();

	for(int k=0;k<N;k++){
		for(int n=0;n<N;n++){
			Xr[k] = Xr[k] + x[n]*cos(2*M_PI*(float)k*(float)n/(float)N);
			Xi[k] = Xi[k] + x[n]*sin(2*M_PI*(float)k*(float)n/(float)N);
		}
		X[k] = sqrt(Xr[k] * Xr[k] + Xi[k]*Xi[k]);
		fprintf(fout,"%f\n",X[k]);
	}
	fclose(fout);
	
	return 0;
}