#include <stdio.h>
#define _USE_MATH_DEFINES
#include <math.h>

float x[10000] = {0,};
float s[10000] = {0,};
float S[10000] = {0,};
float w[10000] = {0,};

float c[10000] = {0,};
int COUNT = 0;
void getAudio(){
	FILE	*fi;
	float	input;
	short	data;

	fi = fopen("Male.raw", "rb");
	for( ; ; ) {
		if(fread(&data, 2, 1, fi) == NULL)
			break;

		input = (float)data;  // DFT에서는 array로 저장
		x[COUNT] = input;
		COUNT++;
		if(COUNT==320)
			break;
	}
	fclose(fi);
}
void setWindow(){
	for(int i=0;i<320;i++){
		w[i]=0.54-0.46*cos(2*M_PI/319*i);
	}
}
void setS(){
	for(int i=0;i<320;i++){
		s[i] = x[i] * w[i];
	}
}
void dft(){
	float Xr[10000] = {0,};
	float Xi[10000] = {0,};

	int N = 320;

	FILE* fout;
	fout = fopen("lpc.txt", "w");

	for(int k=0;k<N;k++){
		for(int n=0;n<N;n++){
			Xr[k] = Xr[k] + s[n]*cos(2*M_PI*(float)k*(float)n/(float)N);
			Xi[k] = Xi[k] + s[n]*sin(2*M_PI*(float)k*(float)n/(float)N);
		}
		S[k] = sqrt(Xr[k] * Xr[k] + Xi[k]*Xi[k]);
		fprintf(fout,"%f\n",S[k]);
	}
	fclose(fout);
}
void dftOfCn(){
	float Xr[10000] = {0,};
	float Xi[10000] = {0,};

	int N = 320;

	FILE* fout;
	fout = fopen("spectralEnvelop.txt", "w");

	for(int k=0;k<N;k++){
		for(int n=0;n<N;n++){
			Xr[k] = Xr[k] + c[n]*cos(2*M_PI*(float)k*(float)n/(float)N);
			Xi[k] = Xi[k] + c[n]*sin(2*M_PI*(float)k*(float)n/(float)N);
		}
		S[k] = sqrt(Xr[k] * Xr[k] + Xi[k]*Xi[k]);
		fprintf(fout,"%f\n",S[k]);
	}
	fclose(fout);
}
void inverseDft(){
	float Xr[10000] = {0,};
	float Xi[10000] = {0,};

	int N = 320;

	FILE* fout;
	fout = fopen("cn.txt", "w");

	for(int k=0;k<N;k++){
		for(int n=0;n<N;n++){
			Xr[k] = Xr[k] + S[n]*cos(2*M_PI*(float)k*(float)n/(float)N);
			Xi[k] = Xi[k] + S[n]*sin(2*M_PI*(float)k*(float)n/(float)N);
		}
		//Xr[k] = Xr[k] / N;
		//Xi[k] = Xi[k] / N;
		c[k] = sqrt(Xr[k] * Xr[k] + Xi[k]*Xi[k]);
		fprintf(fout,"%f\n",c[k]);
	}
	fclose(fout);
}
void liftering(){
	FILE* fout;
	fout = fopen("liftering.txt", "w");
	//liftering
	for(int i=0;i<320;i++){
		if(i > 15 && i < 306){
			c[i] = 0;
		}
		fprintf(fout,"%f\n",c[i]);
	}
	fclose(fout);
}
int main(){
	getAudio();
	setWindow();
	setS();
	dft();

	FILE* fout;
	fout = fopen("cepstrum.txt", "w");
	//S를 log|S|로 변환
	for(int i=0;i<320;i++){
		S[i] = log10(S[i]);
		fprintf(fout,"%f\n",S[i]);
	}

	inverseDft();//cn
	liftering();//vn
	dftOfCn();
	fclose(fout);
	return 0;
}