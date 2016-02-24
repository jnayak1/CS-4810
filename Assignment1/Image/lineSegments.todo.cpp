#include "lineSegments.h"
#include <math.h>

float OrientedLineSegment::length(void) const
{
	float xDisplacement = this->x2 - this->x1;
	float yDisplacement = this->y2 - this->y1;
	float length = (float) sqrt(xDisplacement*xDisplacement + yDisplacement*yDisplacement);

	return length;
}
float OrientedLineSegment::distance(const int& x,const int& y) const
{
	float xDisplacement = this->x2 - this->x1;
	float yDisplacement = this->y2 - this->y1;

	float distance = fabs(yDisplacement*x - xDisplacement*y + (this->x2 * this->y1) + (this->y2 * this->x1)) / this->length();

	return distance;
}
void  OrientedLineSegment::getPerpendicular(float& x,float &y) const
{
	float xDisplacement = this->x2 - this->x1;
	float yDisplacement = this->y2 - this->y1;

	x = (yDisplacement + this->x1) / (this->length());
	y = ((-1*xDisplacement) + this->y1) / (this->length());
}

void  OrientedLineSegment::GetSourcePosition(const OrientedLineSegment& source,const OrientedLineSegment& destination,
											 const int& targetX,const int& targetY,
											 float& sourceX,float& sourceY)
{
	float u, v;
	float uxDiff, uyDiff, vxDiff, vyDiff;

	float xDestDisplacement = destination.x2 - destination.x1;
	float yDestDisplacement = destination.y2 - destination.y1;

	float xDestPerpDisplacementNorm = ((-1) * yDestDisplacement) / destination.length();
	float yDestPerpDisplacementNorm = xDestDisplacement / destination.length();

	u = (((targetX - destination.x1) * (xDestDisplacement)) + ((targetY - destination.y1) * yDestDisplacement)) / (destination.length() * destination.length());
	v = ((targetX - destination.x1) * xDestPerpDisplacementNorm) + ((targetY - destination.y1) * yDestPerpDisplacementNorm);

	uxDiff = u * xDestDisplacement;
	uyDiff = u * yDestDisplacement;

	vxDiff = xDestPerpDisplacementNorm * v;
	vyDiff = yDestPerpDisplacementNorm * v;

	sourceX = destination.x1 + uxDiff + vxDiff;
	sourceY = destination.y1 + uyDiff + vyDiff;
}