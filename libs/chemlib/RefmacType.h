#ifndef REFMAC_TYPE_H
#define REFMAC_TYPE_H

namespace RefmacType
{
    const int MAXLENGTH = 4;        //# of characters in a Refmac type

    const unsigned int DUM  = 0;
    const unsigned int CSP  = 1;    //Sp hybridized, no hydrogen
    const unsigned int CSP1 = 2;    //Triple-bonded, with one hydrogen
    const unsigned int C    = 3;    //Sp2 hybridized, no hydrogen
    const unsigned int C1   = 4;    //Sp2 hybridized, 1 hydrogen
    const unsigned int C2   = 5;    //Sp2 hybridized, 2 hydrogens
    const unsigned int CR1  = 6;
    const unsigned int CR2  = 7;
    const unsigned int CR1H = 8;
    const unsigned int CR15 = 9;    //Sp2 hybridized, in 5-membered ring, 1 hydrogen
    const unsigned int CR5  = 10;   //Sp2 hybridized, in 5-membered ring, no hydrogens
    const unsigned int CR56 = 11;   //Sp2 hybridized, at ring fusion of five- and six-membered rings
    const unsigned int CR55 = 12;   //Sp2 hybridized, at ring fusion of 2 five-membered rings
    const unsigned int CR16 = 13;   //Sp2 hybridized, in 6-membered ring, with one hydrogen
    const unsigned int CR6  = 14;   //Sp2 hybridized, in 6-membered ring, with no hydrogens
    const unsigned int CR66 = 15;   //Sp2 hybridized, at ring fusion of 2 six-membered rings
    const unsigned int CH1  = 16;   //Sp3 hybridized, 1 hydrogen
    const unsigned int CH2  = 17;   //Sp3 hybridized, 2 hydrogens
    const unsigned int CH3  = 18;   //Sp3 hybridized, 3 hydrogens
    const unsigned int CT   = 19;   //Sp3 hybridized, no hydrogens
    const unsigned int NS   = 20;   //Sp hybridized, no hydrogen
    const unsigned int N    = 21;   //Sp2 hybridized, nohydrogen
    const unsigned int NC1  = 22;
    const unsigned int NH1  = 23;
    const unsigned int NC2  = 24;
    const unsigned int NH2  = 25;
    const unsigned int NC3  = 26;
    const unsigned int NT   = 27;
    const unsigned int NT1  = 28;
    const unsigned int NT2  = 29;
    const unsigned int NT3  = 30;
    const unsigned int NPA  = 31;
    const unsigned int NPB  = 32;
    const unsigned int NR5  = 33;
    const unsigned int NR15 = 34;
    const unsigned int NRD5 = 35;
    const unsigned int NR56 = 36;
    const unsigned int NR55 = 37;
    const unsigned int NR6  = 38;
    const unsigned int NR66 = 39;
    const unsigned int NR16 = 40;
    const unsigned int NRD6 = 41;
    const unsigned int OS   = 42;
    const unsigned int O    = 43;
    const unsigned int O2   = 44;
    const unsigned int OH1  = 45;
    const unsigned int OH2  = 46;
    const unsigned int OHA  = 47;
    const unsigned int OHB  = 48;
    const unsigned int OHC  = 49;
    const unsigned int OC2  = 50;
    const unsigned int OC   = 51;
    const unsigned int OP   = 52;
    const unsigned int OB   = 53;
    const unsigned int P    = 54;
    const unsigned int P1   = 55;
    const unsigned int PS   = 56;
    const unsigned int S    = 57;
    const unsigned int S3   = 58;
    const unsigned int S2   = 59;
    const unsigned int S1   = 60;
    const unsigned int ST   = 61;
    const unsigned int SH1  = 62;
    const unsigned int H    = 63;
    const unsigned int HCH  = 64;
    const unsigned int HCH1 = 65;
    const unsigned int HCH2 = 66;
    const unsigned int HCH3 = 67;
    const unsigned int HCR1 = 68;
    const unsigned int HC1  = 69;
    const unsigned int HC2  = 70;
    const unsigned int HCR5 = 71;
    const unsigned int HCR6 = 72;
    const unsigned int HNC1 = 73;
    const unsigned int HNC2 = 74;
    const unsigned int HNH1 = 75;
    const unsigned int HNH2 = 76;
    const unsigned int HNR5 = 77;
    const unsigned int HNR6 = 78;
    const unsigned int HNT1 = 79;
    const unsigned int HNT2 = 80;
    const unsigned int HNT3 = 81;
    const unsigned int HOH1 = 82;
    const unsigned int HOH2 = 83;
    const unsigned int HOHA = 84;
    const unsigned int HOHB = 85;
    const unsigned int HOHC = 86;
    const unsigned int HSH1 = 87;
    const unsigned int SI   = 88;
    const unsigned int SI1  = 89;
    const unsigned int GE   = 90;
    const unsigned int GE1  = 91;
    const unsigned int SN   = 92;
    const unsigned int PB   = 93;
    const unsigned int LI   = 94;
    const unsigned int NA   = 95;
    const unsigned int K    = 96;
    const unsigned int RB   = 97;
    const unsigned int CS   = 98;
    const unsigned int BE   = 99;
    const unsigned int MG   = 100;
    const unsigned int CA   = 101;
    const unsigned int SR   = 102;
    const unsigned int BA   = 103;
    const unsigned int SC   = 104;
    const unsigned int Y    = 105;
    const unsigned int LA   = 106;
    const unsigned int CE   = 107;
    const unsigned int PR   = 108;
    const unsigned int ND   = 109;
    const unsigned int SM   = 110;
    const unsigned int EU   = 111;
    const unsigned int GD   = 112;
    const unsigned int TB   = 113;
    const unsigned int DY   = 114;
    const unsigned int HO   = 115;
    const unsigned int ER   = 116;
    const unsigned int TM   = 117;
    const unsigned int YB   = 118;
    const unsigned int LU   = 119;
    const unsigned int AC   = 120;
    const unsigned int TH   = 121;
    const unsigned int PA   = 122;
    const unsigned int U    = 123;
    const unsigned int NP   = 124;
    const unsigned int TI   = 125;
    const unsigned int ZR   = 126;
    const unsigned int HF   = 127;
    const unsigned int V    = 128;
    const unsigned int NB   = 129;
    const unsigned int TA   = 130;
    const unsigned int CR   = 131;
    const unsigned int MO   = 132;
    const unsigned int W    = 133;
    const unsigned int MN   = 134;
    const unsigned int TC   = 135;
    const unsigned int RE   = 136;
    const unsigned int FE   = 137;
    const unsigned int RU   = 138;
    const unsigned int OSE  = 139;
    const unsigned int CO   = 140;
    const unsigned int RH   = 141;
    const unsigned int IR   = 142;
    const unsigned int NI   = 143;
    const unsigned int PD   = 144;
    const unsigned int PT   = 145;
    const unsigned int CU   = 146;
    const unsigned int AG   = 147;
    const unsigned int AU   = 148;
    const unsigned int ZN   = 149;
    const unsigned int CD   = 150;
    const unsigned int HG   = 151;
    const unsigned int B    = 152;
    const unsigned int AL   = 153;
    const unsigned int GA   = 154;
    const unsigned int IN   = 155;
    const unsigned int TL   = 156;
    const unsigned int AS   = 157;
    const unsigned int AS1  = 158;
    const unsigned int SB   = 159;
    const unsigned int BI   = 160;
    const unsigned int SE   = 161;
    const unsigned int TE   = 162;
    const unsigned int PO   = 163;
    const unsigned int F    = 164;
    const unsigned int CL   = 165;
    const unsigned int BR   = 166;
    const unsigned int I    = 167;
    const unsigned int AT   = 168;
    const unsigned int HE   = 169;
    const unsigned int NE   = 170;
    const unsigned int AR   = 171;
    const unsigned int KR   = 172;
    const unsigned int XE   = 173;
    const unsigned int RN   = 174;
    const unsigned int MAXTYPE = 175;
}

#endif //REFMAC_TYPE_H
