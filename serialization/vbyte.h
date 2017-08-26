#ifndef vbyte_h
#define vbyteh
// Encodes an integer using variable length byte code
// The string where the vbyte is written should be long enough 
// to hold the encoded integers. The function returns the 
// length of the encoding
inline unsigned int vbyte_encode(unsigned long long x,unsigned char * output_string)
{
	unsigned int i;
	for(i=0;;i++)
	{
		output_string[i]=0x80|(x&0x7f);
		x>>=7;
		if(x==0)
		{
			output_string[i]&=0x7f;
			break;

		};
	};
//	printf("Encode length is %d\n",i+1);
	return i+1;
};
// Decodes an integers previously encoded using variable length code

inline unsigned int vbyte_decode(unsigned long long * _x,unsigned char * input_string)
{
	unsigned int i;
	unsigned long long x=0;
	for(i=0;;i++)
	{
		x|=(input_string[i]&0x7f)<<(i*7);
		if((input_string[i]&0x80)==0 || i==9)
			break;
	};
	(*_x)=x;
//	printf("Decode length is %d\n",i+1);
	return i+1;
};


#endif
