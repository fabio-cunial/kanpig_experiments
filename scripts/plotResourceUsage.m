SUMMARY_FILE='../results/summary.txt';


figure(1);


# ------------------------------- FASTQ INPUT ----------------------------------
subplot(1,2,1);
hold on;

A=load(SUMMARY_FILE);
[nRows,nColumns]=size(A);

# Including PBMM2 into sniffles, cutesv, kanpig.
A(2,1)=A(2,1)+A(1,1);
A(3,1)=A(3,1)+A(1,1);
A(6,1)=A(6,1)+A(1,1);
A(2,2)=max([A(2,2),A(1,2)]);
A(3,2)=max([A(3,2),A(1,2)]);
A(6,2)=max([A(6,2),A(1,2)]);

# Plotting
for i=[1:nRows]
    plot(A(i,2),A(i,1),'o');
end
legend('pbmm2','sniffles','cutesv','jedi','minigraph','kanpig','location','northeastoutside');
xlabel('Max RSS (kb)'); ylabel('Wall clock time (s)');
title('FASTQ input');
axis square; grid on;
set(gca,'fontsize',14);



# -------------------------------- BAM INPUT -----------------------------------
subplot(1,2,2);
hold on;

A=load(SUMMARY_FILE);
[nRows,nColumns]=size(A);

plot(A(2,2),A(2,1),'o');
plot(A(3,2),A(3,1),'o');
plot(A(6,2),A(6,1),'o');
legend('sniffles','cutesv','kanpig','location','northeast');
xlabel('Max RSS (kb)'); ylabel('Wall clock time (s)');
title('BAM input');
axis square; grid on;
set(gca,'fontsize',14);
