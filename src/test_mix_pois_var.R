w=c(0.1,0.2,0.3,0.4)
y1=rpois(10000,1)
y2=rpois(10000,7)
y3=rpois(10000,14)
y4=rpois(10000,20)

z1=rpois(10000,5)
z2=rpois(10000,10)
z3=rpois(10000,15)
z4=rpois(10000,20)


x1=w[1]*y1+w[2]*y2+w[3]*y3+w[4]*y4
x2=w[1]*z1+w[2]*z2+w[3]*z3+w[4]*z4



mean(x1/(x1+x2))
var(x1/(x1+x2))

m1=w[1]*1+w[2]*7+w[3]*14+w[4]*20
m2=w[1]*5+w[2]*10+w[3]*15+w[4]*20

m1/(m1+m2)
mean(((m1/(m1+m2))*(m2/(m1+m2))/(x1+x2)))



y=rpois(10000,1)
z=rpois(10000,24)


mean(y/(y+z),na.rm=TRUE)
var(y/(y+z),na.rm=TRUE)


1/25
mean(1/25*24/25/(y+z)[(y+z)!=0])
