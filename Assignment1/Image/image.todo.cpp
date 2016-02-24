#include "image.h"
#include <stdlib.h>
#include <math.h>
#include <cmath>
#include <climits>
#include <iostream>

// allows the use of min() and max() functions
#include <algorithm>
using namespace std;


Pixel::Pixel(const Pixel32& p)
{
}

Pixel32::Pixel32(const Pixel& p)
{
}

int Image32::AddRandomNoise(const float& noise,Image32& outputImage) const
{
	int height = this->height();
	int width = this->width();
	outputImage.setSize(width, height);
	int noiseRange = (noise * 256) * 2;
	
	for(int y=0; y<height; ++y){
		for(int x=0; x<width; ++x){
			int randNum;
			if(noiseRange != 0) randNum = (rand() % noiseRange) - (.5*noiseRange);
			else{ randNum = 0; }
			outputImage.pixel(x,y).r = max(0,min(255,randNum + this->pixel(x,y).r));
			outputImage.pixel(x,y).g = max(0,min(255,randNum + this->pixel(x,y).g));
			outputImage.pixel(x,y).b = max(0,min(255,randNum + this->pixel(x,y).b));
			outputImage.pixel(x,y).a = this->pixel(x,y).a;
		}	
	}
	return 1;	
}

int Image32::Brighten(const float& brightness,Image32& outputImage) const
{
	int height = this->height();
	int width = this->width();
	outputImage.setSize(width, height);

	for(int y=0; y<height; ++y){
		for(int x=0; x<width; ++x){
			int r = this->pixel(x,y).r;
			int g = this->pixel(x,y).g;
			int b = this->pixel(x,y).b;
			int rBright = int(r * brightness);
			int gBright = int(g * brightness);
			int bBright = int(b * brightness);
			if(rBright <= 255 && gBright <= 255 && bBright <= 255 ){
				outputImage.pixel(x,y).r = rBright;
				outputImage.pixel(x,y).g = gBright;
				outputImage.pixel(x,y).b = bBright;
				outputImage.pixel(x,y).a = this->pixel(x,y).a;
			}
			else{
				// must clamp if any r,g,b values are over 255
				int maxBrightness = 255 / max(r,max(g,b));
				outputImage.pixel(x,y).r = int(r * maxBrightness);
				outputImage.pixel(x,y).g = int(g * maxBrightness);
				outputImage.pixel(x,y).b = int(b * maxBrightness);
				outputImage.pixel(x,y).a = this->pixel(x,y).a;
			}
		}
	}
	return 1;
}

int Image32::Luminance(Image32& outputImage) const
{
	int height = this->height();
	int width = this->width();
	outputImage.setSize(width, height);

	for(int y=0; y<height; ++y){
		for(int x=0; x<width; ++x){
			int r = this->pixel(x,y).r;
			int g = this->pixel(x,y).g;
			int b = this->pixel(x,y).b;

			int l = int(0.30*r + 0.59*g + 0.11*b);

			outputImage.pixel(x,y).r = l;
			outputImage.pixel(x,y).g = l;
			outputImage.pixel(x,y).b = l;
			outputImage.pixel(x,y).a = this->pixel(x,y).a;
		}
	}	
	return 1;
}

int Image32::Contrast(const float& contrast,Image32& outputImage) const
{
	int height = this->height();
	int width = this->width();
	outputImage.setSize(width, height);

	float totalLum = 0.0;
	int numPixels = 0;

	for(int y=0; y<height; ++y){
		for(int x=0; x<width; ++x){
			int r = this->pixel(x,y).r;
			int g = this->pixel(x,y).g;
			int b = this->pixel(x,y).b;
			float l = 0.30*r + 0.59*g + 0.11*b;
			totalLum += l;
			numPixels += 1;
		}
	}	
	float meanLum = totalLum / numPixels;

	for(int y=0; y<height; ++y){
		for(int x=0; x<width; ++x){
			int r = this->pixel(x,y).r;
			int g = this->pixel(x,y).g;
			int b = this->pixel(x,y).b;
			int a = this->pixel(x,y).a;

			float l = 0.30*r + 0.59*g + 0.11*b;
			float lDiff = l - meanLum;

			float maxContrast;
			if(lDiff >= 0){ 
				maxContrast = (255 - (max(r,max(g,b)) - lDiff)) / lDiff;
			}
			else{
				maxContrast = fabs(((min(r,min(g,b))) - lDiff) / lDiff); // lDiff is negative, so adding to minimum of r,g,b
			} 

			float contrastClamped = min(maxContrast,contrast);

			outputImage.pixel(x,y).r = contrastClamped*lDiff + (r-lDiff);
			outputImage.pixel(x,y).g = contrastClamped*lDiff + (g-lDiff);
			outputImage.pixel(x,y).b = contrastClamped*lDiff + (b-lDiff);
			outputImage.pixel(x,y).a = a;
		}
	}	
	return 1;
}

int Image32::Saturate(const float& saturation,Image32& outputImage) const
{
	int height = this->height();
	int width = this->width();
	outputImage.setSize(width, height);

	for(int y=0; y<height; ++y){
		for(int x=0; x<width; ++x){
			int r = this->pixel(x,y).r;
			int g = this->pixel(x,y).g;
			int b = this->pixel(x,y).b;

			float l = 0.30*r + 0.59*g + 0.11*b;

			float rDiff = r-l;
			float gDiff = g-l;
			float bDiff = b-l;

			bool rOverflow = (r-rDiff) + rDiff*saturation > 255 || (r-rDiff) + rDiff*saturation < 0;
			bool gOverflow = (g-gDiff) + gDiff*saturation > 255 || (g-gDiff) + gDiff*saturation < 0;
			bool bOverflow = (b-bDiff) + bDiff*saturation > 255 || (b-bDiff) + bDiff*saturation < 0;

			float satClamped = saturation; // if no overflows, satClamped is saturation
			if(rOverflow){ // if r overflows, clamp sat to max so r doesn't overflow
				satClamped = min((r-rDiff), 255-(r-rDiff)) / rDiff;
			}
			if(gOverflow){ // if g overflows, clamp sat to max so g doesn't overflow (or satClamped if less)
				satClamped = min(min((g-gDiff), 255-(g-gDiff)) / gDiff, satClamped);
			}
			if(bOverflow){ // if b overflows, clamp sat to max so b doesn't overflow (or satClamped if less)
				satClamped = min(min((b-bDiff), 255-(b-bDiff)) / bDiff, satClamped);
			}
			outputImage.pixel(x,y).r = satClamped*rDiff + (r-rDiff);
			outputImage.pixel(x,y).g = satClamped*gDiff + (g-gDiff);
			outputImage.pixel(x,y).b = satClamped*bDiff + (b-bDiff);
			outputImage.pixel(x,y).a = this->pixel(x,y).a;
		}
	}	
	return 1;
}

int Image32::Quantize(const int& bits,Image32& outputImage) const
{
	int height = this->height();
	int width = this->width();
	outputImage.setSize(width, height);

	for(int y=0; y<height; ++y){
		for(int x=0; x<width; ++x){
			int r = this->pixel(x,y).r;
			int g = this->pixel(x,y).g;
			int b = this->pixel(x,y).b;

			int multiple = int(pow(2,8-bits)); 

			int quanR = r + multiple/2;
			quanR = min(255, quanR - (quanR % multiple)); // clamped to 255

			int quanG = g + multiple/2;
			quanG = min(255, quanG - (quanG % multiple)); // clamped to 255

			int quanB = b + multiple/2;
			quanB = min(255, quanB - (quanB % multiple)); // clamped to 255
			
			outputImage.pixel(x,y).r = quanR;
			outputImage.pixel(x,y).g = quanG;
			outputImage.pixel(x,y).b = quanB;
			outputImage.pixel(x,y).a = this->pixel(x,y).a;
		}
	}
	return 1;
}

int Image32::RandomDither(const int& bits,Image32& outputImage) const
{
	float noise = 0.5/pow(2,bits-1);
	Image32 interImage;
	this->AddRandomNoise(noise, interImage);
	interImage.Quantize(bits,outputImage);
	return 1;
}

int Image32::OrderedDither2X2(const int& bits,Image32& outputImage) const
{
	int height = this->height();
	int width = this->width();
	int i;
	int j;
	float D[2][2];
	D[0][0] = 1.0;
	D[0][1] = 3.0;
	D[1][0] = 4.0;
	D[1][1] = 2.0;

	Image32 interImage;
	interImage.setSize(width,height);

	float scaler255 = (255.0/pow(2,bits-1));

	this->Quantize(bits,interImage);


	for(int y=0; y<height; ++y){
		for (int x = 0; x<width; ++x){
			i = x % 2;
			j = y % 2;

			float r = float(this->pixel(x,y).r) / 255.0;
			float g = float(this->pixel(x,y).g) / 255.0;
			float b = float(this->pixel(x,y).b) / 255.0;

			float cR = r*pow(2,bits-1);
			float cG = g*pow(2,bits-1);
			float cB = b*pow(2,bits-1);

			float eR =  cR - floor(cR);
			float eG =  cG - floor(cG);
			float eB =  cB - floor(cB);

			if(eR > (D[i][j] / 8)) {outputImage.pixel(x,y).r = int(ceil(cR)*scaler255);}
			else{outputImage.pixel(x,y).r = int(floor(cR) * scaler255);}

			if(eG > (D[i][j] / 8)) {outputImage.pixel(x,y).g = int(ceil(cG)*scaler255);}
			else{outputImage.pixel(x,y).g = int(floor(cG) * scaler255);}

			if(eB > (D[i][j] / 8)) {outputImage.pixel(x,y).b = int(ceil(cB)*scaler255);}
			else{outputImage.pixel(x,y).b = int(floor(cB) * scaler255);}

			outputImage.pixel(x,y).a = this->pixel(x,y).a;
		}
	}
	return 1;
}

int Image32::FloydSteinbergDither(const int& bits,Image32& outputImage) const
{
	int height = this->height();
	int width = this->width();
	float alpha = 7/16;
	float beta = 3/16;
	float gamma = 5/16;
	float delta = 1/16;

	this->Quantize(bits,outputImage);
	for (int y = 0; y < height; ++y){
		for (int x = 0; x < width; ++x){
			float eR = this->pixel(x,y).r - outputImage.pixel(x,y).r;
			float eG = this->pixel(x,y).g - outputImage.pixel(x,y).g;
			float eB = this->pixel(x,y).b - outputImage.pixel(x,y).b;

			if(y+1 < height){
				outputImage.pixel(x,y+1).r = max(0,min(255,int(outputImage.pixel(x,y+1).r + alpha * eR)));
				outputImage.pixel(x,y+1).g = max(0,min(255,int(outputImage.pixel(x,y+1).g + alpha * eG)));
				outputImage.pixel(x,y+1).b = max(0,min(255,int(outputImage.pixel(x,y+1).b + alpha * eB)));
			}
			if(x+1 < width && y > 0){
				outputImage.pixel(x+1,y-1).r = max(0,min(255,int(outputImage.pixel(x+1,y-1).r + beta * eR)));
				outputImage.pixel(x+1,y-1).g = max(0,min(255,int(outputImage.pixel(x+1,y-1).g + beta * eG)));
				outputImage.pixel(x+1,y-1).b = max(0,min(255,int(outputImage.pixel(x+1,y-1).b + beta * eB)));
			}
			if(x+1 < width){
				outputImage.pixel(x+1,y).r = max(0,min(255,int(outputImage.pixel(x+1,y).r + gamma * eR)));
				outputImage.pixel(x+1,y).g = max(0,min(255,int(outputImage.pixel(x+1,y).g + gamma * eG)));
				outputImage.pixel(x+1,y).b = max(0,min(255,int(outputImage.pixel(x+1,y).b + gamma * eB)));
			}
			if(x+1 < width && y+1 < height){
				outputImage.pixel(x+1,y+1).r = max(0,min(255,int(outputImage.pixel(x+1,y+1).r + delta * eR)));
				outputImage.pixel(x+1,y+1).g = max(0,min(255,int(outputImage.pixel(x+1,y+1).g + delta * eG)));
				outputImage.pixel(x+1,y+1).b = max(0,min(255,int(outputImage.pixel(x+1,y+1).b + delta * eB)));
			}
		}
	}
	return 1;
}

int Image32::Blur3X3(Image32& outputImage) const
{
	int width = this->width();
	int height = this->height();

	for(int y=0; y<height; ++y){
		for(int x=0; x<width; ++x){
			float denominator;

			bool leftBorder = (x==0);
			bool rightBorder = (x == width-1);
			bool bottomBorder = (y == height-1);
			bool xBorder = leftBorder || rightBorder;
			bool yBorder = ( y==0 || y == (height-1));

			if(xBorder || yBorder){				
				outputImage.pixel(x,y).r = this->pixel(x,y).r;
				outputImage.pixel(x,y).g = this->pixel(x,y).g;
				outputImage.pixel(x,y).b = this->pixel(x,y).b;
			}
			else{
				denominator = 16;

				float topR = (this->pixel(x-1,y-1).r + 2*this->pixel(x,y-1).r  + this->pixel(x+1,y-1).r) / (denominator);
				float middleR = (2*this->pixel(x-1,y).r + 4*this->pixel(x,y).r  + 2*this->pixel(x+1,y).r) / (denominator);
				float bottomR = (this->pixel(x-1,y+1).r + 2*this->pixel(x,y+1).r  + this->pixel(x+1,y+1).r) / (denominator);
				outputImage.pixel(x,y).r = topR + middleR + bottomR;

				float topG = (this->pixel(x-1,y-1).g + 2*this->pixel(x,y-1).g  + this->pixel(x+1,y-1).g) / (denominator);
				float middleG = (2*this->pixel(x-1,y).g + 4*this->pixel(x,y).g  + 2*this->pixel(x+1,y).g) / (denominator);
				float bottomG = (this->pixel(x-1,y+1).g + 2*this->pixel(x,y+1).g  + this->pixel(x+1,y+1).g) / (denominator);
				outputImage.pixel(x,y).g = topG + middleG + bottomG;

				float topB = (this->pixel(x-1,y-1).b + 2*this->pixel(x,y-1).b  + this->pixel(x+1,y-1).b) / (denominator);
				float middleB = (2*this->pixel(x-1,y).b + 4*this->pixel(x,y).b  + 2*this->pixel(x+1,y).b) / (denominator);
				float bottomB = (this->pixel(x-1,y+1).b + 2*this->pixel(x,y+1).b  + this->pixel(x+1,y+1).b) / (denominator);
				outputImage.pixel(x,y).b = topB + middleB + bottomB;
			}
			outputImage.pixel(x,y).a = this->pixel(x,y).a;
		}
	}


	return 1;
}

int Image32::EdgeDetect3X3(Image32& outputImage) const
{
	int width = this->width();
	int height = this->height();

	for(int y=0; y<height; ++y){
		for(int x=0; x<width; ++x){
			float denominator;

			bool leftBorder = (x==0);
			bool rightBorder = (x == width-1);
			bool bottomBorder = (y == height-1);
			bool xBorder = leftBorder || rightBorder;
			bool yBorder = ( y==0 || y == (height-1));

			if(xBorder || yBorder){				
				outputImage.pixel(x,y).r = this->pixel(x,y).r;
				outputImage.pixel(x,y).g = this->pixel(x,y).g;
				outputImage.pixel(x,y).b = this->pixel(x,y).b;
			}
			else{
				denominator = 16;

				float topR = (this->pixel(x-1,y-1).r + this->pixel(x,y-1).r  + this->pixel(x+1,y-1).r) / (-1*denominator);
				float middleR = (-1*this->pixel(x-1,y).r + 8*this->pixel(x,y).r  + -1*this->pixel(x+1,y).r) / (denominator);
				float bottomR = (this->pixel(x-1,y+1).r + this->pixel(x,y+1).r  + this->pixel(x+1,y+1).r) / (-1*denominator);
				outputImage.pixel(x,y).r = topR + middleR + bottomR;

				float topG = (this->pixel(x-1,y-1).g + this->pixel(x,y-1).g  + this->pixel(x+1,y-1).g) / (-1*denominator);
				float middleG = (-1*this->pixel(x-1,y).g + 8*this->pixel(x,y).g  + -1*this->pixel(x+1,y).g) / (denominator);
				float bottomG = (this->pixel(x-1,y+1).g + this->pixel(x,y+1).g  + this->pixel(x+1,y+1).g) / (-1*denominator);
				outputImage.pixel(x,y).g = topG + middleG + bottomG;

				float topB = (this->pixel(x-1,y-1).b + this->pixel(x,y-1).b  + this->pixel(x+1,y-1).b) / (-1*denominator);
				float middleB = (-1*this->pixel(x-1,y).b + 8*this->pixel(x,y).b  + -1*this->pixel(x+1,y).b) / (denominator);
				float bottomB = (this->pixel(x-1,y+1).b + this->pixel(x,y+1).b  + this->pixel(x+1,y+1).b) / (-1*denominator);
				outputImage.pixel(x,y).b = topB + middleB + bottomB;
			}
			outputImage.pixel(x,y).a = this->pixel(x,y).a;
		}
	}
	return 1;
}
int Image32::ScaleNearest(const float& scaleFactor,Image32& outputImage) const
{
	int widthSRC = this->width();
	int heightSRC = this->height();

	int widthDST = int(floor(widthSRC * scaleFactor + 0.5));
	int heightDST = int(floor(heightSRC * scaleFactor + 0.5));

	outputImage.setSize(widthDST, heightDST);

	for(int yDST=0; yDST<heightDST-(int)scaleFactor; ++yDST){
		for(int xDST=0; xDST<widthDST-(int)scaleFactor; ++xDST){
			outputImage.pixel(xDST, yDST) = this->NearestSample(((float)xDST / scaleFactor), ((float)yDST / scaleFactor));
		}
	}
	return 1;
}

int Image32::ScaleBilinear(const float& scaleFactor,Image32& outputImage) const
{
	int widthSRC = this->width();
	int heightSRC = this->height();

	int widthDST = int(floor(widthSRC * scaleFactor + 0.5));
	int heightDST = int(floor(heightSRC * scaleFactor + 0.5));

	outputImage.setSize(widthDST, heightDST);

	for(int yDST=0; yDST<heightDST-(int)scaleFactor; ++yDST){
		for(int xDST=0; xDST<widthDST-(int)scaleFactor; ++xDST){
			outputImage.pixel(xDST, yDST) = this->BilinearSample(((float)xDST / scaleFactor),((float)yDST / scaleFactor));
		}
	}
	return 1;
}

int Image32::ScaleGaussian(const float& scaleFactor,Image32& outputImage) const
{
	int widthSRC = this->width();
	int heightSRC = this->height();

	int widthDST = int(floor((float)widthSRC * scaleFactor)-ceil(scaleFactor));
	int heightDST = int(floor((float)heightSRC * scaleFactor)-ceil(scaleFactor));


	outputImage.setSize(widthDST, heightDST);

	
	float variance = 1.0 / scaleFactor;
	float r = 3.0;

	for(int yDST=0; yDST < heightDST; ++yDST){
		for(int xDST=0; xDST < widthDST; ++xDST){
			outputImage.pixel(xDST, yDST) = this->GaussianSample(((float)xDST / scaleFactor), ((float)yDST / scaleFactor), variance, r);
		}
	}
	return 1;
}

int Image32::RotateNearest(const float& angle,Image32& outputImage) const
{
	int height = this->height();
	int width = this->width();

	double rad = -1.0 * angle * (M_PI/180.0);

	double midHeight = (height/2.0);
	double midWidth = (width/2.0);

	double diagonal = sqrt(midHeight*midHeight + midWidth*midWidth);

	int offset = ceil(diagonal - min(midWidth, midHeight));

	Image32 interImage;
	interImage.setSize(width+(2*offset), height+(2*offset));

	outputImage.setSize(width+(2*offset), height+(2*offset));

	

	for(int y=0; y < height; ++y){
		for(int x=0; x < width; ++x){
			int xOffset = x + offset;
			int yOffset = y + offset;

			interImage.pixel(xOffset, yOffset).r = this->pixel(x,y).r;
			interImage.pixel(xOffset, yOffset).g = this->pixel(x,y).g;
			interImage.pixel(xOffset, yOffset).b = this->pixel(x,y).b;
		}
	}

	for(int y=-1*offset; y < height + offset; ++y){
		for(int x=-1*offset; x < width + offset; ++x){
			int xOffset = x + offset;
			int yOffset = y + offset;

			float u = cos(rad)*(x - midWidth) - sin(rad)*(y - midHeight) + offset + midWidth;
			float v = sin(rad)*(x - midWidth) + cos(rad)*(y - midHeight) + offset + midHeight;


			Pixel32 pix = interImage.NearestSample(u,v);


			outputImage.pixel(xOffset,yOffset).r = pix.r;
			outputImage.pixel(xOffset,yOffset).g = pix.g;
			outputImage.pixel(xOffset,yOffset).b = pix.b;
			outputImage.pixel(xOffset,yOffset).a = pix.a;

			// useful for debugging
			// outputImage.pixel(xOffset,yOffset).r = interImage.pixel(xOffset, yOffset).r;
			// outputImage.pixel(xOffset,yOffset).g = interImage.pixel(xOffset, yOffset).g;
			// outputImage.pixel(xOffset,yOffset).b = interImage.pixel(xOffset, yOffset).b;
		}
	}

	return 1;
}

int Image32::RotateBilinear(const float& angle,Image32& outputImage) const
{
	int height = this->height();
	int width = this->width();

	double rad = -1.0 * angle * (M_PI/180.0);

	double midHeight = (height/2.0);
	double midWidth = (width/2.0);

	double diagonal = sqrt(midHeight*midHeight + midWidth*midWidth);

	int offset = ceil(diagonal - min(midWidth, midHeight));

	Image32 interImage;
	interImage.setSize(width+(2*offset), height+(2*offset));

	outputImage.setSize(width+(2*offset), height+(2*offset));

	

	for(int y=0; y < height; ++y){
		for(int x=0; x < width; ++x){
			int xOffset = x + offset;
			int yOffset = y + offset;

			interImage.pixel(xOffset, yOffset).r = this->pixel(x,y).r;
			interImage.pixel(xOffset, yOffset).g = this->pixel(x,y).g;
			interImage.pixel(xOffset, yOffset).b = this->pixel(x,y).b;
		}
	}

	for(int y= (-1*offset); y < height + offset; ++y){
		for(int x= (-1*offset); x < width + offset; ++x){
			int xOffset = x + offset;
			int yOffset = y + offset;

			cout << rad << " " << x << " " << y << "\n";

			float u = cos(rad)*(x - midWidth) - sin(rad)*(y - midHeight) + offset + midWidth;
			float v = sin(rad)*(x - midWidth) + cos(rad)*(y - midHeight) + offset + midHeight;

			cout << u << " " << v << "\n";

			Pixel32 pix = interImage.BilinearSample(u,v);


			outputImage.pixel(xOffset,yOffset).r = pix.r;
			outputImage.pixel(xOffset,yOffset).g = pix.g;
			outputImage.pixel(xOffset,yOffset).b = pix.b;
			outputImage.pixel(xOffset,yOffset).a = pix.a;


			// useful for debugging
			// outputImage.pixel(xOffset,yOffset).r = interImage.pixel(xOffset, yOffset).r;
			// outputImage.pixel(xOffset,yOffset).g = interImage.pixel(xOffset, yOffset).g;
			// outputImage.pixel(xOffset,yOffset).b = interImage.pixel(xOffset, yOffset).b;
		}
	}

	return 1;
}
	
int Image32::RotateGaussian(const float& angle,Image32& outputImage) const
{
	int height = this->height();
	int width = this->width();

	double rad = -1.0 * angle * (M_PI/180.0);

	double midHeight = (height/2.0);
	double midWidth = (width/2.0);

	double diagonal = sqrt(midHeight*midHeight + midWidth*midWidth);

	int offset = ceil(diagonal - min(midWidth, midHeight));

	Image32 interImage;
	interImage.setSize(width+(2*offset), height+(2*offset));

	outputImage.setSize(width+(2*offset), height+(2*offset));

	float variance = 2;
	float radius = 2;

	for(int y=0; y < height; ++y){
		for(int x=0; x < width; ++x){
			int xOffset = x + offset;
			int yOffset = y + offset;

			interImage.pixel(xOffset, yOffset).r = this->pixel(x,y).r;
			interImage.pixel(xOffset, yOffset).g = this->pixel(x,y).g;
			interImage.pixel(xOffset, yOffset).b = this->pixel(x,y).b;
		}
	}

	for(int y=-1*offset; y < height + offset; ++y){
		for(int x=-1*offset; x < width + offset; ++x){
			int xOffset = x + offset;
			int yOffset = y + offset;

			float u = cos(rad)*(x - midWidth) - sin(rad)*(y - midHeight) + offset + midWidth;
			float v = sin(rad)*(x - midWidth) + cos(rad)*(y - midHeight) + offset + midHeight;


			Pixel32 pix = interImage.GaussianSample(u,v, variance, radius);


			outputImage.pixel(xOffset,yOffset).r = pix.r;
			outputImage.pixel(xOffset,yOffset).g = pix.g;
			outputImage.pixel(xOffset,yOffset).b = pix.b;
			outputImage.pixel(xOffset,yOffset).a = pix.a;

			// useful for debugging
			// outputImage.pixel(xOffset,yOffset).r = interImage.pixel(xOffset, yOffset).r;
			// outputImage.pixel(xOffset,yOffset).g = interImage.pixel(xOffset, yOffset).g;
			// outputImage.pixel(xOffset,yOffset).b = interImage.pixel(xOffset, yOffset).b;
		}
	}

	return 1;
}


int Image32::SetAlpha(const Image32& matte)
{
	int height = this->height();
	int width = this->width();
	for(int y=0; y<height; ++y){
		for(int x=0; x<width; ++x){
			this->pixel(x,y).a = matte.pixel(x,y).r; // one channel, all components should be equal
		}
	}
	return 1;
}

int Image32::Composite(const Image32& overlay,Image32& outputImage) const
{
	int height = this->height();
	int width = this->width();
	int heightOverlay = overlay.height();
	int widthOverlay = overlay.width();

	outputImage.setSize(width, height);

	for(int y=0; y<height; ++y){
		for(int x=0; x<width; ++x){
			int alpha256 = overlay.pixel(x,y).a;
			float alpha = (float) alpha256 / 255.0;
			outputImage.pixel(x,y).r = alpha*(overlay.pixel(x,y).r) + (1.0-alpha)*(this->pixel(x,y).r);
			outputImage.pixel(x,y).g = alpha*(overlay.pixel(x,y).g) + (1.0-alpha)*(this->pixel(x,y).g);
			outputImage.pixel(x,y).b = alpha*(overlay.pixel(x,y).b) + (1.0-alpha)*(this->pixel(x,y).b);
			outputImage.pixel(x,y).a = 255;
		}
	}
	return 1;
}

int Image32::CrossDissolve(const Image32& source,const Image32& destination,const float& blendWeight,Image32& outputImage)
{
	int height = destination.height();
	int width = destination.width();


	cout << "CrossDissolve" << "\n";
	for(int y=0; y<height; ++y){
		for(int x=0; x<width; ++x){
			outputImage.pixel(x,y).r = (blendWeight * source.pixel(x,y).r) +  ((1-blendWeight) * destination.pixel(x,y).r);
			outputImage.pixel(x,y).g = (blendWeight * source.pixel(x,y).g) +  ((1-blendWeight) * destination.pixel(x,y).g);
			outputImage.pixel(x,y).b = (blendWeight * source.pixel(x,y).b) +  ((1-blendWeight) * destination.pixel(x,y).b);
		}
	}
	
	return 1;
}
int Image32::Warp(const OrientedLineSegmentPairs& olsp,Image32& outputImage) const
{
	int height = this->height();
	int width = this->width();
	outputImage.setSize(width,height);

	int numOfLineSegments = olsp.count;

	float dSumX, dSumY, weight, weightSum;
	cout << "warp" << "\n";
	for(int y=0; y<height; ++y){
		for(int x=0; x<width; ++x){
			for(int i=0; i<numOfLineSegments; ++i){
				dSumX = 0, dSumY = 0, weightSum = 0;
				weight = olsp.segments2[i].getWeight(x,y);
				weightSum += weight;

				float sourceX, sourceY;

				olsp.segments1[i].GetSourcePosition(olsp.segments1[i], olsp.segments2[i], x, y, sourceX, sourceY);

				dSumX += (sourceX - x) * weight;
				dSumY += (sourceY - y) * weight;
			}
			
			 
			outputImage.pixel(x,y).r = this->pixel(x + (dSumX/weightSum), y + (dSumY/weightSum)).r;
			outputImage.pixel(x,y).g = this->pixel(x + (dSumX/weightSum), y + (dSumY/weightSum)).g;
			outputImage.pixel(x,y).b = this->pixel(x + (dSumX/weightSum), y + (dSumY/weightSum)).b;
		}
	}
	return 1;
}

int Image32::FunFilter(Image32& outputImage) const
{
	double height = this->height();
	double width = this->width();

	for(int y=0; y<height; ++y){
		for(int x=0; x<width; ++x){
				double y1 =  (double) y /  (double) (height-1.0);
				double x1 = (double) x / (double) (width-1.0);
				double r = sqrt( pow(x1 - 0.5, 2.0) + pow(y1 - 0.5, 2.0) );

				double a = atan2((y1 - 0.5) , (x1 - 0.5));
				
				double rd = pow(r, 2.5)/.5;

				float xNew = ((rd*cos(a)) + 0.5) * (width-1);
				float yNew = ((rd*sin(a)) + 0.5) * (height-1);

				Pixel32 pix = this->NearestSample(xNew,yNew);

				outputImage.pixel(x,y).r = pix.r;
				outputImage.pixel(x,y).g = pix.g;
				outputImage.pixel(x,y).b = pix.b;
		}
	}

	return 1;
}

int Image32::Crop(const int& x1,const int& y1,const int& x2,const int& y2,Image32& outputImage) const
{
	int width = (x2 - x1) + 1;
	int height = (y2 - y1) + 1;
	outputImage.setSize(width,height);
	int i = 0;
	int j = 0;

	for(int y=y1; y <= y2; ++y){
		i = 0;
		for(int x=x1; x <= x2; ++x){
			outputImage.pixel(i,j).r = this->pixel(x,y).r;
			outputImage.pixel(i,j).g = this->pixel(x,y).g;
			outputImage.pixel(i,j).b = this->pixel(x,y).b;
			outputImage.pixel(i,j).a = this->pixel(x,y).a;
			i += 1;
		}
		j += 1;
	}
	return 1;
}

Pixel32 Image32::NearestSample(const float& x,const float& y) const
{
	int height = this->height();
	int width = this->width();
	Pixel32 nearestPixel;
	int nearestX = min(width-1,max(0,int(floor(x+0.5))));
	int nearestY = min(height-1,max(0,int(floor(y+0.5))));
	nearestPixel.r = this->pixel(nearestX,nearestY).r;
	nearestPixel.g = this->pixel(nearestX,nearestY).g;
	nearestPixel.b = this->pixel(nearestX,nearestY).b;
	nearestPixel.a = this->pixel(nearestX,nearestY).a;

	return nearestPixel;
}

Pixel32 Image32::BilinearSample(const float& x,const float& y) const
{
	Pixel32 bilinearPixel;
	int height = this->height();
	int width = this->width();


	int left = min(width-1, max(0,int(floor(x))));
	int right = min(width-1, max(0,int(ceil(x))));
	int up =  min(height-1, max(0,int(floor(y))));
	int down =  min(height-1, max(0,int(ceil(y))));

	float dx = x - floor(x);
	float dy = y - floor(y);


	
	float aR = this->pixel(left,down).r * (1-dx) + this->pixel(right,down).r * dx;
	float bR = this->pixel(left,up).r * (1-dx) + this->pixel(right,up).r * dx;

	float aG = this->pixel(left,down).g * (1-dx) + this->pixel(right,down).g * dx;
	float bG = this->pixel(left,up).g * (1-dx) + this->pixel(right,up).g * dx;

	float aB = this->pixel(left,down).b * (1-dx) + this->pixel(right,down).b * dx;
	float bB = this->pixel(left,up).b * (1-dx) + this->pixel(right,up).b * dx;

	bilinearPixel.r = aR*(1-dy) + bR*dy;
	bilinearPixel.g = aG*(1-dy) + bG*dy;
	bilinearPixel.b = aB*(1-dy) + bB*dy;

	return bilinearPixel;
}

Pixel32 Image32::GaussianSample(const float& x,const float& y,const float& variance,const float& radius) const
{
	int width = this->width();
	int height = this->height();

	Pixel32 gaussianPixel;

	float gaussTotalR = 0.0;
	float gaussTotalG = 0.0;
	float gaussTotalB = 0.0;

	float denominator = 0.0;

	for(int j=1; j<=radius; ++j){
		for(int i=1; i<=radius; ++i){
			double left = ceil(x - i);
			double right = floor(x + i);
			double top = ceil(x - i);
			double bottom = floor(x + i);

			cout << left << " " << right << " " << top << " " << bottom;

			double topLeftDist = sqrt((left-x)*(left-x) + (top-y)*(top-y));
			double topRightDist = sqrt((right-x)*(right-x) + (top-y)*(top-y));
			double bottomLeftDist = sqrt((left-x)*(left-x) + (bottom-y)*(bottom-y));
			double bottomRightDist = sqrt((right-x)*(right-x) + (bottom-y)*(bottom-y));

			double gaussTopLeft = exp(-1.0*(left*left + top*top)/(2.0*variance*variance));
			double gaussTopRight = exp(-1.0*(right*right + top*top)/(2.0*variance*variance));
			double gaussBottomLeft = exp(-1.0*(left*left + bottom*bottom)/(2.0*variance*variance));
			double gaussBottomRight = exp(-1.0*(right*right + bottom*bottom)/(2.0*variance*variance));

			if( (topLeftDist < radius) && (top >= 0) && (left >= 0) && (top < height) && (left < width)){
				gaussTotalR += gaussTopLeft * (float)this->pixel(left,top).r; 
				gaussTotalG += gaussTopLeft * (float)this->pixel(left,top).g; 
				gaussTotalB += gaussTopLeft * (float)this->pixel(left,top).b;

				cout << "hi";

				denominator += gaussTopLeft;
			}
			if( (topRightDist < radius) && (top >= 0) && (right < width) && (bottom < top) && right >= 0 ){
				gaussTotalR += gaussTopRight * (float)this->pixel(right,top).r; 
				gaussTotalG += gaussTopRight * (float)this->pixel(right,top).g;
				gaussTotalB += gaussTopRight * (float)this->pixel(right,top).b;

				denominator += gaussTopRight;
			}
			if( (bottomLeftDist < radius) && (bottom < height) && (left >=0) && left < width && bottom >= 0){
				gaussTotalR += gaussBottomLeft * (float)this->pixel(left,bottom).r;
				gaussTotalG += gaussBottomLeft * (float)this->pixel(left,bottom).g;
				gaussTotalB += gaussBottomLeft * (float)this->pixel(left,bottom).b;

				denominator += gaussBottomLeft;
			}
			if( (bottomRightDist < radius) && (bottom < height) && (right < width) && bottom >= 0 && right >=0){
				gaussTotalR += gaussBottomRight * (float)this->pixel(right,bottom).r;
				gaussTotalG += gaussBottomRight * (float)this->pixel(right,bottom).g;
				gaussTotalB += gaussBottomRight * (float)this->pixel(right,bottom).b;
				denominator += gaussBottomRight;
			}
			cout<< "\n ";
		}

		

		gaussianPixel.r = int(floor((gaussTotalR / denominator) + 0.5));
		gaussianPixel.g = int(floor((gaussTotalG / denominator) + 0.5));
		gaussianPixel.b = int(floor((gaussTotalB / denominator) + 0.5));
	}

	return gaussianPixel;
}
