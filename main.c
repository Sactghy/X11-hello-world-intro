#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <X11/Xlib.h>
#include <X11/Xutil.h>

const int dts = 4122+6, x_c = 1020, y_c = 740, lnn = 390, lnm = 660, lncol = 0, lncoV = 0, lncoB = 0;

void circle(char *dt, double rad, int x, int y, unsigned char dens, unsigned char mask)
{	unsigned char tr, bp,gp,rp;
    double xc = 0.0, xs, yc = 1.0, rc = rad*rad;
    int pos = 0, k1 = 1, k2 = 1;
    int cc = 1, cr = 1;
    while ( yc >= 0 ) { k2 = 1; if ( cr ) { xs = xc*xc; yc = rc-xs; yc = sqrt(yc); }
     else { xs = yc*yc; xc = rc-xs; xc = sqrt(xc); }
     if ( ( cr == 1 ) && ( (int)xc >= (int)yc ) ) { cr = 0; yc = xc;}
     for ( int b = 1; b <= 4; b++) { if ( ( xc > x ) || ( xc < 0 ) ) goto sk1;
      pos = ( ( ( x_c * ( y + ( k1 * (int)yc ) ) ) + ( x + ( k2 * (int)xc ) ) ) * 4 );
      if ( (pos <= 0 ) || ( pos >= x_c * y_c * 4 ) ) goto sk1; if ( dt[pos+3] == mask ) goto sk1;

       if (x_c==2000) { unsigned char ma=dt[pos+3];

        tr=dt[pos+2]; if ( dens > tr || (int)(tr-dens)<0)tr=dens; else tr &= dens;

        if ( ( ma & 0b10000000 ) && !( ma & 0b01000000 ) ) {
        if (lncol==0) { bp=tr-tr/4; gp=tr-tr/2; rp=tr-tr/3; }
        if (lncol==1) { bp=tr-tr/4; gp=tr; rp=tr-tr/2; }
        if (lncol==2) { bp=tr-tr/2; gp=tr-tr/3; rp=tr; } }

        if ( ( ma & 0b01000000 ) && !( ma & 0b10000000 ) ) {
        if (lncoV==0) { bp=tr-tr/2; gp=tr-tr/3; rp=tr-tr; }
        if (lncoV==1) { bp=1; gp=tr-tr/5; rp=tr+tr/5; }
        if (lncoV==2) { bp=tr; gp=tr-tr; rp=tr; } }

        tr = dt[pos+0]; if ( dens > tr || (int)(tr+dens) < 0 ) tr = dens; else tr &= dens; dt[pos+0] = bp;
        tr = dt[pos+1]; if ( dens > tr || (int)(tr+dens) < 0 ) tr = dens; else tr &= dens; dt[pos+1] = gp;
        tr = dt[pos+2]; if ( dens > tr || (int)(tr+dens) < 0 ) tr = dens; else tr &= dens; dt[pos+2] = rp;
        dt[pos+3] = mask; goto sk1;} // end

      tr=dt[pos]; if ( dens > tr ) tr=0; else tr -= dens;
      dt[pos]=tr; dt[pos+1]=tr; dt[pos+2]=tr; dt[pos+3]=mask;

sk1:  if ( ( b % 2 != 0 ) ) k2 = -1; else k2 = 1;
      switch ( b ) { case 2: k1 = -1; break; case 4: k1 = 1; break; }
     } if ( cr ) xc++; else yc--; }
}

void centerc( char *dt, unsigned char dens )
{	dens+=13;
    unsigned char tr;
    int x = 510-2, y = 370-2, pos = 0, diag = 2, flat = 3;
    int xy = 1, yx = 1, ab = 1, cc = 1;
    double xc,xs,yc,rc,r = 20.0;
    for (int i=1; i<=25; i++) { pos=(((1020*y)+x)*4);
      tr=dt[pos]; if (dens>tr) tr=0; else tr -= dens;
      dt[pos]=tr; dt[pos+1]=tr; dt[pos+2]=tr; dt[pos+3]=11; x++;
      if (i%5==0) {y++; x -=5;} } x++; for (int i=1; i<=12; i++) {
      pos=(((1020*y)+x)*4);
      tr=dt[pos]; if (dens>tr) tr=0; else tr -= dens;
      dt[pos]=tr; dt[pos+1]=tr; dt[pos+2]=tr; dt[pos+3]=11;
      (i<=6)?x++:y++; if((i%3==0)&&(i<6)){y-=6;x-=3;} if((i%6==0)){y+=2; x++;}
      if((i%9==0)){x-=6;y-=3;}}
    x++;
    for(int m=0; m<=16; m++){ for(int i=0; i<=3; i++){ for (int n=1; n<=flat; n++)
        { pos=(((1020*y)+x)*4); tr=dt[pos]; if (dens>tr) tr=0; else tr -= dens;
      dt[pos]=tr; dt[pos+1]=tr; dt[pos+2]=tr; dt[pos+3]=11;
      if(i<=1)(xy)?x++:y--;else(xy)?x--:y++; } if(i<=1)(xy)?y--:x--;else(xy)?y++:x++;
      for(int n=1; n<=diag; n++){ pos=(((1020*y)+x)*4);
      tr=dt[pos]; if (dens>tr) tr=0; else tr -= dens;
      dt[pos]=tr; dt[pos+1]=tr; dt[pos+2]=tr; dt[pos+3]=11;
          if(i<=1){if(xy){x++; y--;} else {x--; y--;}}else{if(xy){x--; y++;} else {x++; y++;}}}
      xy = !xy; } y++; diag++;
      if (diag>flat) { y--; x += flat;
      for(int i=0; i<=3; i++){ for(int n=1; n<=diag; n++){ pos=(((1020*y)+x)*4);
      tr=dt[pos]; if (dens>tr) tr=0; else tr -= dens;
      dt[pos]=tr; dt[pos+1]=tr; dt[pos+2]=tr; dt[pos+3]=11;
      if(i==0){x++; y--;} if(i==1){x--; y--;} if(i==2){x--; y++;} if(i==3){x++; y++;}}
      if(i==0){x--; y-=flat;} if(i==1){y++; x-=flat;} if(i==2){x++; y+=flat;} }
      x--; flat+=2; diag--; dens -= 2; }
    } dens += 2;
nw:     circle(dt,r,510,370,dens,11); r+=0.5; dens-=1; if (dens>3) goto nw; }

struct Circo
{
    int slx;
    double c1_rad;
    unsigned char c2_grow;
};

void initC(struct Circo *c)
{ c->c1_rad = 48; c->slx = 18; c->c2_grow = 19; }

void CircGo( char *dt, struct Circo *c )
{	unsigned char gr1, gr0; double rr1 = c->c1_rad,rr2 = c->c1_rad;
    circle(dt, c->c1_rad, 510, 370, c->c2_grow, 64);
    gr1 = c->c2_grow;
ag_:	gr1 -= 1; if (gr1>c->c2_grow) goto sk_;
    gr0 = gr1-c->slx; if (gr0>c->c2_grow) goto sk_;
    rr1+=0.5; rr2-=0.5;
    circle(dt, rr1, 510, 370, gr0, 64); circle(dt, rr2, 510, 370, gr0, 64);
    goto ag_;
sk_:	if ((c->slx!=0)&&(c->c2_grow%18==0)) c->slx -=3;
    if (c->c1_rad<715) c->c1_rad++; else {c->c1_rad=44;c->c2_grow=19;c->slx=18; goto nx_;}
    if (c->c1_rad>30)   c->c1_rad +=0.4; if (c->c1_rad>60)   c->c1_rad +=0.5;
    if (c->c1_rad>80)  {c->c1_rad +=0.5; if ((c->c2_grow<235)) {c->c2_grow++;}}
    if (c->c1_rad>95)  {c->c1_rad +=0.6; if ((c->c2_grow<235)) {c->c2_grow++;}}
    if (c->c1_rad>135)  c->c1_rad +=0.7; if (c->c1_rad>165)  c->c1_rad +=0.9;
    if (c->c1_rad>195)  c->c1_rad +=1.1; if (c->c1_rad>225)  c->c1_rad +=1.3;
    if (c->c1_rad>255)  c->c1_rad +=1.9;
nx_:	if (c->c2_grow<144) c->c2_grow++;
}

struct Randm{ unsigned char res, result; };

char rand0( struct Randm* a )
{	a->res = rand();
    a->result = (unsigned char)(a->res);
    return a->result; 	 }

int main(int argc, char *argv[])
{
    srand(423); rand();
    time_t now;
    struct tm beg = *localtime(&now);
    double seconds;

    Display *d = XOpenDisplay(0);

    XSetWindowAttributes swa;
    swa.event_mask = KeyPressMask|ButtonPressMask|ButtonReleaseMask|ButtonMotionMask;
    if ( d )
    { Window w = XCreateWindow(d, XDefaultRootWindow(d), 0, 0, 1018, 740, 0, 			CopyFromParent, CopyFromParent, CopyFromParent, 0b00100000000000, &swa);
    XStoreName(d, w, "Hello world\0");
    GC gc = XDefaultGC(d, 0);
    XMapWindow(d, w); XFlush (d);
//----------------------------------------------------------------
    int x_r = 0, y_r = 0, st1, st2, st3, rWaX, rWaY;
    unsigned int w_r = 0, h_r = 0, b_w_r = 0, d_r = 0, nc_r;
    Window root_r, parent_r, child_r;
    XWindowAttributes rWa;

bad_st:	st1 = XGetWindowAttributes(d, w, &rWa);
    printf(" %i %i : rWa1\n", rWa.x, rWa.y );

    Window *ccW = &child_r;
    st2 = XQueryTree(d, w, &root_r, &parent_r, &ccW, &nc_r);
    XGetGeometry(d, parent_r, &child_r, &x_r, &y_r, &w_r, &h_r, &b_w_r, &d_r);
    printf("%i %i %i %i %i %i : Geometry\n",x_r,y_r,w_r,h_r,b_w_r,d_r);

    st3 = XGetWindowAttributes(d, w, &rWa);
    printf("%i %i - %i : rWa2\n", rWa.x, rWa.y, rWa.override_redirect );

    if ((st1==0)||(st2==0)||(st3==0))
    {printf("Bad rWa -> getting position again...\n"); goto bad_st;}

    root_r = XDefaultRootWindow(d);
    int s_x= 0, s_y= 0, dx_r, dy_r;
    if (XTranslateCoordinates(d, root_r, w, s_x, s_y, &dx_r, &dy_r, &child_r))
    { printf("%i %i -xTrC", dx_r, dy_r );
    if (( (x_r+rWa.x) != ((-1)*dx_r) ) && ( (y_r+rWa.y) != ((-1)*dy_r) ) ||
       ( (rWa. x== ((-1)*dx_r) ) && ( rWa.y == ((-1)*dx_r) ))) goto bad_st;} else goto bad_st;

    XSizeHints hnt; hnt.flags = PMaxSize|PMinSize;
    hnt.min_width=1020; hnt.min_height=740; hnt.max_width=1020; hnt.max_height=740;
    XSetWMNormalHints(d, w, &hnt);

    XImage *Img = XGetImage( d, root_r, x_r+rWa.x, y_r+rWa.y, 1020, 740, AllPlanes, ZPixmap );
    if (!Img) { printf("Can't get: Img\n"); }

uu:	; char *d2 = malloc( sizeof(char[1020*740*4]));
    if (!d2) {printf("Bad alloc: d2\n"); goto uu;}
    char *dt = Img->data;
    for (int i=0; i<1020*740*4; i+= 4){ d2[i]=dt[i];d2[i+1]=dt[i+1];d2[i+2]=dt[i+2]; }

        XPutImage( d, w, gc, Img, 0, 0, 0, 0, 1020, 740 );

oo:	; char *dtX = malloc( sizeof(char[2000*2000*4]));
    if (!dtX) {printf( "Bad alloc: dtX\n"); goto oo;}
    XImage *Img2 = XCreateImage(d, XDefaultVisual(d,0), Img->depth, ZPixmap, Img->xoffset, dtX, 2000, 2000, Img->bitmap_pad, 8000);
    if ( !Img2 ) {printf("Can't create: Img2 \n"); }
    for (int i=0; i<=2000*2000*4; i += 4) { dtX[i]=0;dtX[i+1]=0;dtX[i+2]=0;}

o3:	; char *dtM = malloc( sizeof (char[1020*740*4]));
    if (!dtM) {printf( "Bad alloc: dtM\n" ); goto o3;}
    XImage *Img4 = XCreateImage(d, XDefaultVisual(d,0), Img->depth, ZPixmap, Img->xoffset, dtM, 1020, 740, Img->bitmap_pad, 4080);
    if (!Img4) { printf("Can't create: Img4\n"); }

    XEvent e, e1; int x = 1, average = 0, seg = 1, pos = 0, cor= 0;
    int stage1= 0, st2chk = 0, bpress= 0, stage3 = 0, stage4 = 0;
    int back = 0, done= 0, nextstep= 0, yn = 0, is = 0, cond;
    int crcmode = 0, stage5 = 0, scene3 = 0, stars = 0;

    unsigned char r,g,b,sh= 0,shf = 0b11111111, rnd= 0, cnt_dots= 0, col0, cl1= 0;
    unsigned char r1nd, div_r, trns, prec = 32, c_grow = 2;
    int uuc, urc, rrc, drc, ddc, dlc, llc, ulc;
    int cnt0= 0, cnt1= 0, cdd = 0, circcnt= 0, ycnt = 13, xcnt = 17, ycnt2 = 17, xcnt2 = 13;
    struct Circo C1,C2,C3,C4,C5,C6;
    initC(&C1); initC(&C2); initC(&C3); initC(&C4); initC(&C5); initC(&C6);
    struct Randm vRandom;

    char cntcnt3= 0, prcnt3= 0;

    unsigned int px, py, ln1x = 400, chckLn, chckLn2, chckLn3;
    unsigned int* ox = malloc ( sizeof(unsigned int[dts]) );
    unsigned int* oy = malloc ( sizeof(unsigned int[dts]) );
    int* ob = malloc( sizeof(int[dts]));

    do { for (int i=0; i<1020*740*4; i += 4 )
     { if (!crcmode) rnd = rand0(&vRandom); else rnd=123; rnd &= sh;
         if ( st2chk ) {dt[i] = 128+rnd; dt[i+1] = 128+rnd; dt[i+2] = 128+rnd;
                 d2[i] = 128+rnd; goto ex2;}
         if ( (x % 1017000 == 0 ) && ( sh < 127 ) ) {sh++;} if ( sh == 127 ) stage1=1;
       r=(unsigned char)(d2[i+0]);
       g=(unsigned char)(d2[i+1]);
       b=(unsigned char)(d2[i+2]);
        if (r>128) d2[i+0]--; if (r<128) d2[i+0]++;
        if (g>128) d2[i+1]--; if (g<128) d2[i+1]++;
        if (b>128) d2[i+2]--; if (b<128) d2[i+2]++;
        if ((d2[i]+rnd)>255||(d2[i+1]+rnd)>255||(d2[i+2]+rnd)>255) {rnd |= 0b00000111;}
        dt[i] = d2[i]+rnd; dt[i+1] = d2[i+1]+rnd; dt[i+2] = d2[i+2]+rnd;
    ex2: x++;
      }
        if ( st2chk ) for ( int i = 0; i < 1020*740*4; i += 4 ) {
        if ( prec == 0 ) { average=dt[i]; goto skip; }
            average=dt[i]-prec;


    skip:	if ((seg>=2)&&((i<(1020*(187*(seg-1))))*4)||(i>(1020*(740-(46*(seg-1))))*4))
        goto blur1; else {r1nd = d2[i]; trns=0; goto set;}
    blur1:	rrc=i+4; if (rrc<0||rrc>1020*740*4) average +=34; else average +=d2[rrc];
        if ((seg>=3)&&((i<(1020*(184*(seg-2))))*4)||(i>(1020*(740-(46*(seg-2))))*4))
        goto blur2; else {cor=seg-2; goto blend;}
    blur2:	drc=i+(1021*4);
        if (drc<0||drc>1020*740*4) average +=34; else average +=d2[drc];
        if ((seg>=4)&&((i<(1020*(184*(seg-3))))*4)||(i>(1020*(740-(46*(seg-3))))*4))
        goto blur3; else {cor=seg-3; goto blend;}
    blur3:	ddc=i+(1020*4);
        if (ddc<0||ddc>1020*740*4) average +=34; else average +=d2[ddc];
        if ((seg>=5)&&((i<(1020*(184*(seg-4))))*4)||(i>(1020*(740-(46*(seg-4))))*4))
        goto blur4; else {cor=seg-4; goto blend;}
    blur4:	llc=i-4;
        if (llc<0||llc>1020*740*4) average +=34; else average +=d2[llc];
        if ((seg>=6)&&((i<(1020*(184*(seg-5))))*4)||(i>(1020*(740-(46*(seg-5))))*4))
        goto blur5; else {cor=seg-5; goto blend;}
    blur5:	dlc=i+(1019*4);
        if (dlc<0||dlc>1020*740*4) average +=34; else average +=d2[dlc];
        if ((seg>=7)&&((i<(1020*(184*(seg-6))))*4)||(i>(1020*(740-(46*(seg-6))))*4))
        goto blur6; else {cor=seg-6; goto blend;}
    blur6:	uuc=i-(1020*4);
        if (uuc<0||uuc>1020*740*4) average +=34; else average +=d2[uuc];
        if ((seg>=8)&&((i<(1020*(184*(seg-7))))*4)||(i>(1020*(740-(46*(seg-7))))*4))
        goto blur7; else {cor=seg-7; goto blend;}
    blur7:	urc=i-(1019*4);
        if (urc<0||urc>1020*740*4) average +=34; else average +=d2[urc];
        if ((seg>=9)&&((i<(1020*(186*(seg-8))))*4)||(i>(1020*(740-(46*(seg-8))))*4))
        goto blur8; else {cor=seg-8; goto blend;}
    blur8:	ulc=i-(1021*4);
        if (ulc<0||ulc>1020*740*4) average +=34; else average +=d2[ulc];
        cor=seg-9;

    blend:	div_r = average % (seg-cor); average /= (seg-cor); average += (div_r*5);
        trns = (char)(average); r1nd = dt[i]; trns &= ~shf;
        if (!prec) {trns += prec+div_r; r1nd -= prec+div_r;} r1nd &= shf;
    set:
        if ((cnt0%2464800==0)&&(shf!=0b00000000)) {shf<<=1;prec-=4;cdd=cnt0;}
        if ((cnt0==cdd)&&(shf==0b00000000)&&(seg<9)) {seg++;cdd=cnt0+3464800;
                         if (seg==3) stage3=1; if (seg==7) stage4=1;}

        dt[i]=r1nd+trns; dt[i+1]=r1nd+trns; dt[i+2]=r1nd+trns; dt[i+3]=3; cnt0++;
    }
    if (stage1) st2chk=1;
    if (stage3) { if (back) { c_grow -=2; if (c_grow<169) {back = !back;} goto cc_; }
            if (c_grow<255-15) {c_grow +=4; if (c_grow>=255-15) back = !back; }
        cc_:	 centerc(dt, c_grow ); }

    if (stage4) {	if (circcnt!=1000) circcnt++;
                    if (circcnt>36)  CircGo(dt,&C1);
                    if (circcnt>78)  CircGo(dt,&C2);
                    if (circcnt>124) CircGo(dt,&C3);
                    if (circcnt>200) CircGo(dt,&C4);
                    if (circcnt>300) CircGo(dt,&C5);
                    if (circcnt>400) CircGo(dt,&C6);	 }

    int x,y,x_m,y_m;
    if(XCheckWindowEvent(d,w,ButtonReleaseMask|ButtonPressMask|ButtonMotionMask,&e))
    {
      if (e.xbutton.type==ButtonPress) { x=e.xbutton.x; y=e.xbutton.y; bpress=1; x_m=x; y_m=y; }
      if (e.xbutton.type==ButtonRelease) bpress=0;
      if (e.xmotion.type==MotionNotify) { x_m=e.xmotion.x; y_m=e.xmotion.y; }
    }

    if (bpress==1) if ((x-x_m>4)||(x-x_m<-4)&&(y-y_m>4)||(y-y_m<-4))
    { int dif_x, dif_y, kx, ky; unsigned char reverse;
      if (x_m>1019) x_m=1019; if (x_m<0) x_m=0; if (y_m>739) y_m=739; if (y_m<0) y_m=0;
      if ( x>x_m ) {dif_x = x-x_m; kx = -1;} else { dif_x = x_m-x; kx = 1;}
      if ( y>y_m ) {dif_y = y-y_m; ky = -1;} else { dif_y = y_m-y; ky = 1;}
        for(int s=0; s<=dif_x; s++)
          {dt[(((y*1020)+x+(kx*s))*4)+1]=0; dt[(((y*1020)+x+(kx*s))*4)+2]=0;
           dt[((((y+(ky*dif_y))*1020)+x+(kx*s))*4)+1]=0;
            dt[((((y+(ky*dif_y))*1020)+x+(kx*s))*4)+2]=0;}
        for(int s=0; s<=dif_y; s++)
          {dt[((((y+(s*ky))*1020)+x)*4)+1]=0;dt[((((y+(s*ky))*1020)+x)*4)+2]=0;
        dt[((((y+(s*ky))*1020)+(x+(kx*dif_x)))*4)+1]=0;
            dt[((((y+(s*ky))*1020)+(x+(kx*dif_x)))*4)+2]=0;}
        for(int s=0; s<dif_y-1; s++) { for(int u=0; u<dif_x-1; u++)
        { reverse=dt[(((((ky*(s+1))+y)*1020)+(kx*(u+1))+x)*4)+1];
            dt[(((((ky*(s+1))+y)*1020)+(kx*(u+1))+x)*4)+1]=~reverse;
        reverse=dt[(((((ky*(s+1))+y)*1020)+(kx*(u+1))+x)*4)+2];
            dt[(((((ky*(s+1))+y)*1020)+(kx*(u+1))+x)*4)+2]=~reverse; } }

    } XPutImage(d, w, gc, Img, 0, 0, 0, 0, 1020, 740);

    XCheckWindowEvent(d,w,KeyPressMask,&e1);

    if ((e1.type==KeyPress)&&(e1.xkey.keycode==65)) {e1.type=0; XNextEvent(d, &e1); e1.type=0; }

    if (!done&&(e1.type==KeyPress)&&(e1.xkey.keycode==52)) {e1.type=0; crcmode=!crcmode; } }

    while ((e1.type!=KeyPress)&&((e1.xkey.keycode>255)||(e1.xkey.keycode==65)||(e1.xkey.keycode==52)));

    printf( "Key code: %i\n", e1.xkey.keycode);

    free(d2); free(dtX); free(dtM);
    free(ox); free(oy); free(ob);

    } XCloseDisplay(d);
}
