function test_DTW(template,data)

d1=template;
d2=data;

[y1 pla1]=PLA(d1,1,1);
[y2 pla2]=PLA(d2,1,1);
[w ta tb] = simmx_dtw(y1,pla1,y2,pla2);
[p,q,Dm] = dp_dtw(w);
[ym1 ym2 yout1]=draw_dtw(y1,pla1,p,y2,pla2,q);

draw_dtw
