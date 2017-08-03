args<-commandArgs()


#print (args[6]);
#print (args[7]);
#print (args[8]);

x<-as.real(args[6]);
y<-as.real(args[7]);
p<-as.double(args[8]);

result<-"";

result<-binom.test(x,y,p=p,alternative="two.sided")

print(result$p.value);
