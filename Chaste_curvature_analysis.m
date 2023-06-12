close all
clc
clear



B1_0 = readmatrix("../Chaste_connectivity/CC_B1_small_0.csv");
B1_1 = readmatrix("../Chaste_connectivity/CC_B1_small_1.csv");
B1_2 = readmatrix("../Chaste_connectivity/CC_B1_small_2.csv");
B1_3 = readmatrix("../Chaste_connectivity/CC_B1_small_3.csv");
B1_4 = readmatrix("../Chaste_connectivity/CC_B1_small_4.csv");

B1_boundary = [B1_0(:,1);B1_1(:,1);B1_2(:,1);B1_3(:,1);B1_4(:,1) ];
B1_nonboundary = [B1_0(:,3);B1_1(:,3);B1_2(:,3);B1_3(:,3);B1_4(:,3) ];
B1_curvature = [B1_0(:,2);B1_1(:,2);B1_2(:,2);B1_3(:,2);B1_4(:,2) ];



B2_0 = readmatrix("../Chaste_connectivity/CC_B2_small_0.csv");
B2_1 = readmatrix("../Chaste_connectivity/CC_B2_small_1.csv");
B2_2 = readmatrix("../Chaste_connectivity/CC_B2_small_2.csv");
B2_3 = readmatrix("../Chaste_connectivity/CC_B2_small_3.csv");
B2_4 = readmatrix("../Chaste_connectivity/CC_B2_small_4.csv");

B2_boundary = [B2_0(:,1);B2_1(:,1);B2_2(:,1);B2_3(:,1);B2_4(:,1) ];
B2_nonboundary = [B2_0(:,3);B2_1(:,3);B2_2(:,3);B2_3(:,3);B2_4(:,3) ];
B2_curvature = [B2_0(:,2);B2_1(:,2);B2_2(:,2);B2_3(:,2);B2_4(:,2) ];


B3_0 = readmatrix("../Chaste_connectivity/CC_B3_small_0.csv");
B3_1 = readmatrix("../Chaste_connectivity/CC_B3_small_1.csv");
B3_2 = readmatrix("../Chaste_connectivity/CC_B3_small_2.csv");
B3_3 = readmatrix("../Chaste_connectivity/CC_B3_small_3.csv");
B3_4 = readmatrix("../Chaste_connectivity/CC_B3_small_4.csv");

B3_boundary = [B3_0(:,1);B3_1(:,1);B3_2(:,1);B3_3(:,1);B3_4(:,1) ];
B3_nonboundary = [B3_0(:,3);B3_1(:,3);B3_2(:,3);B3_3(:,3);B3_4(:,3) ];
B3_curvature = [B3_0(:,2);B3_1(:,2);B3_2(:,2);B3_3(:,2);B3_4(:,2) ];



B4_0 = readmatrix("../Chaste_connectivity/CC_B4_small_0.csv");
B4_1 = readmatrix("../Chaste_connectivity/CC_B4_small_1.csv");
B4_2 = readmatrix("../Chaste_connectivity/CC_B4_small_2.csv");
B4_3 = readmatrix("../Chaste_connectivity/CC_B4_small_3.csv");
B4_4 = readmatrix("../Chaste_connectivity/CC_B4_small_4.csv");

B4_boundary = [B4_0(:,1);B4_1(:,1);B4_2(:,1);B4_3(:,1);B4_4(:,1) ];
B4_nonboundary = [B4_0(:,3);B4_1(:,3);B4_2(:,3);B4_3(:,3);B4_4(:,3) ];
B4_curvature = [B4_0(:,2);B4_1(:,2);B4_2(:,2);B4_3(:,2);B4_4(:,2) ];



figure;
subplot(1,4,1)
scatter(B1_boundary./B1_nonboundary,B1_curvature);
ylabel("$\kappa$")
xlabel("$n_{1}/n_{2}$")
title("B1")
ylim([-0.2,0.2])

subplot(1,4,2)
scatter(B2_boundary./B2_nonboundary,B2_curvature);
ylabel("$\kappa$")
xlabel("$n_{1}/n_{2}$")
title("B2")
ylim([-0.2,0.2])

subplot(1,4,3)
scatter(B3_boundary./B3_nonboundary,B3_curvature);
ylabel("$\kappa$")
xlabel("$n_{1}/n_{2}$")
title("B3")
ylim([-0.2,0.2])

subplot(1,4,4)
scatter(B4_boundary./B4_nonboundary,B4_curvature);
ylabel("$\kappa$")
xlabel("$n_{1}/n_{2}$")
title("B4")
ylim([-0.2,0.2])


%

all_bound = [B1_boundary;B2_boundary;B3_boundary;B4_boundary];
all_nonbound = [B1_nonboundary;B2_nonboundary;B3_nonboundary;B4_nonboundary];
all_curv = [B1_curvature;B2_curvature;B3_curvature;B4_curvature];






figure;
scatter(all_bound ./all_nonbound ,all_curv);
%B1
n2_1_B1_curv = B1_curvature(B1_nonboundary(:) == 1);
n2_2_B1_curv = B1_curvature(B1_nonboundary(:) == 2);
n2_3_B1_curv = B1_curvature(B1_nonboundary(:) == 3);
n2_4_B1_curv = B1_curvature(B1_nonboundary(:) == 4);

n2_1_B1_curv  = [ n2_1_B1_curv; NaN(length(n2_2_B1_curv)-length(n2_1_B1_curv),1)];
n2_2_B1_curv  = [ n2_2_B1_curv; NaN(length(n2_2_B1_curv)-length(n2_2_B1_curv),1)];
n2_3_B1_curv  = [ n2_3_B1_curv; NaN(length(n2_2_B1_curv)-length(n2_3_B1_curv),1)];
n2_4_B1_curv  = [ n2_4_B1_curv; NaN(length(n2_2_B1_curv)-length(n2_4_B1_curv),1)];


%B2
n2_1_B2_curv = B2_curvature(B2_nonboundary(:) == 1);
n2_2_B2_curv = B2_curvature(B2_nonboundary(:) == 2);
n2_3_B2_curv = B2_curvature(B2_nonboundary(:) == 3);
n2_4_B2_curv = B2_curvature(B2_nonboundary(:) == 4);

n2_1_B2_curv  = [ n2_1_B2_curv; NaN(length(n2_2_B2_curv)-length(n2_1_B2_curv),1)];
n2_2_B2_curv  = [ n2_2_B2_curv; NaN(length(n2_2_B2_curv)-length(n2_2_B2_curv),1)];
n2_3_B2_curv  = [ n2_3_B2_curv; NaN(length(n2_2_B2_curv)-length(n2_3_B2_curv),1)];
n2_4_B2_curv  = [ n2_4_B2_curv; NaN(length(n2_2_B2_curv)-length(n2_4_B2_curv),1)];

%B3
n2_1_B3_curv = B3_curvature(B3_nonboundary(:) == 1);
n2_2_B3_curv = B3_curvature(B3_nonboundary(:) == 2);
n2_3_B3_curv = B3_curvature(B3_nonboundary(:) == 3);
n2_4_B3_curv = B3_curvature(B3_nonboundary(:) == 4);

n2_1_B3_curv  = [ n2_1_B3_curv; NaN(length(n2_2_B3_curv)-length(n2_1_B3_curv),1)];
n2_2_B3_curv  = [ n2_2_B3_curv; NaN(length(n2_2_B3_curv)-length(n2_2_B3_curv),1)];
n2_3_B3_curv  = [ n2_3_B3_curv; NaN(length(n2_2_B3_curv)-length(n2_3_B3_curv),1)];
n2_4_B3_curv  = [ n2_4_B3_curv; NaN(length(n2_2_B3_curv)-length(n2_4_B3_curv),1)];

%B4
n2_1_B4_curv = B4_curvature(B4_nonboundary(:) == 1);
n2_2_B4_curv = B4_curvature(B4_nonboundary(:) == 2);
n2_3_B4_curv = B4_curvature(B4_nonboundary(:) == 3);
n2_4_B4_curv = B4_curvature(B4_nonboundary(:) == 4);

n2_1_B4_curv  = [ n2_1_B4_curv; NaN(length(n2_2_B4_curv)-length(n2_1_B4_curv),1)];
n2_2_B4_curv  = [ n2_2_B4_curv; NaN(length(n2_2_B4_curv)-length(n2_2_B4_curv),1)];
n2_3_B4_curv  = [ n2_3_B4_curv; NaN(length(n2_2_B4_curv)-length(n2_3_B4_curv),1)];
n2_4_B4_curv  = [ n2_4_B4_curv; NaN(length(n2_2_B4_curv)-length(n2_4_B4_curv),1)];



figure();
subplot(1,4,1)
violinplot([n2_1_B1_curv,n2_2_B1_curv,n2_3_B1_curv,NaN(length(n2_2_B1_curv),1)],["1","2","3","4"]);
   ylabel("Average Curvature $(\mu m^{-1})$")
   xlabel("Number of heterotypic connections, $n_{2}$")
   title("Homotypic connections $n_{1} = 2$")
ylim([-0.3,0.2])
xlim([0,5])
   
   subplot(1,4,2)
violinplot([n2_1_B2_curv,n2_2_B2_curv,n2_3_B2_curv,n2_4_B2_curv],["1","2","3","4"]);
ylim([-0.3,0.2])
xlim([0,5])

subplot(1,4,3)
violinplot([n2_1_B3_curv,n2_2_B3_curv,n2_3_B3_curv,n2_4_B3_curv],["1","2","3","4"]);
ylim([-0.3,0.2])
xlim([0,5])

subplot(1,4,4)
violinplot([n2_1_B4_curv,n2_2_B4_curv,n2_3_B4_curv,n2_4_B4_curv],["1","2","3","4"]);
ylim([-0.3,0.2])
xlim([0,5])

%% medium size organoid

close all 
clc
clear


B1_0 = readmatrix("../Chaste_connectivity/CC_B1_med_0.csv");
B1_1 = readmatrix("../Chaste_connectivity/CC_B1_med_1.csv");
B1_2 = readmatrix("../Chaste_connectivity/CC_B1_med_2.csv");
B1_3 = readmatrix("../Chaste_connectivity/CC_B1_med_3.csv");
B1_4 = readmatrix("../Chaste_connectivity/CC_B1_med_4.csv");

B1_5 = readmatrix("../Chaste_connectivity/CC_B1_med_5.csv");
B1_6 = readmatrix("../Chaste_connectivity/CC_B1_med_6.csv");
B1_7 = readmatrix("../Chaste_connectivity/CC_B1_med_7.csv");
B1_8 = readmatrix("../Chaste_connectivity/CC_B1_med_8.csv");
B1_9 = readmatrix("../Chaste_connectivity/CC_B1_med_9.csv");

B1_boundary = [B1_0(:,1);B1_1(:,1);B1_2(:,1);B1_3(:,1);B1_4(:,1);
                B1_5(:,1);B1_6(:,1);B1_7(:,1);B1_8(:,1);B1_9(:,1)];
            
B1_nonboundary = [B1_0(:,3);B1_1(:,3);B1_2(:,3);B1_3(:,3);B1_4(:,3);
                B1_5(:,3);B1_6(:,3);B1_7(:,3);B1_8(:,3);B1_9(:,3) ];

B1_curvature = [B1_0(:,2);B1_1(:,2);B1_2(:,2);B1_3(:,2);B1_4(:,2);
                B1_5(:,2);B1_6(:,2);B1_7(:,2);B1_8(:,2);B1_9(:,2) ];



B2_0 = readmatrix("../Chaste_connectivity/CC_B2_med_0.csv");
B2_1 = readmatrix("../Chaste_connectivity/CC_B2_med_1.csv");
B2_2 = readmatrix("../Chaste_connectivity/CC_B2_med_2.csv");
B2_3 = readmatrix("../Chaste_connectivity/CC_B2_med_3.csv");
B2_4 = readmatrix("../Chaste_connectivity/CC_B2_med_4.csv");

B2_5 = readmatrix("../Chaste_connectivity/CC_B2_med_5.csv");
B2_6 = readmatrix("../Chaste_connectivity/CC_B2_med_6.csv");
B2_7 = readmatrix("../Chaste_connectivity/CC_B2_med_7.csv");
B2_8 = readmatrix("../Chaste_connectivity/CC_B2_med_8.csv");
B2_9 = readmatrix("../Chaste_connectivity/CC_B2_med_9.csv");


B2_boundary = [B2_0(:,1);B2_1(:,1);B2_2(:,1);B2_3(:,1);B2_4(:,1);
                B2_5(:,1);B2_6(:,1);B2_7(:,1);B2_8(:,1);B2_9(:,1)];
            
B2_nonboundary = [B2_0(:,3);B2_1(:,3);B2_2(:,3);B2_3(:,3);B2_4(:,3);
                B2_5(:,3);B2_6(:,3);B2_7(:,3);B2_8(:,3);B2_9(:,3) ];
            
B2_curvature = [B2_0(:,2);B2_1(:,2);B2_2(:,2);B2_3(:,2);B2_4(:,2);
                B2_5(:,2);B2_6(:,2);B2_7(:,2);B2_8(:,2);B2_9(:,2) ];


B3_0 = readmatrix("../Chaste_connectivity/CC_B3_med_0.csv");
B3_1 = readmatrix("../Chaste_connectivity/CC_B3_med_1.csv");
B3_2 = readmatrix("../Chaste_connectivity/CC_B3_med_2.csv");
B3_3 = readmatrix("../Chaste_connectivity/CC_B3_med_3.csv");
B3_4 = readmatrix("../Chaste_connectivity/CC_B3_med_4.csv");

B3_5 = readmatrix("../Chaste_connectivity/CC_B3_med_5.csv");
B3_6 = readmatrix("../Chaste_connectivity/CC_B3_med_6.csv");
B3_7 = readmatrix("../Chaste_connectivity/CC_B3_med_7.csv");
B3_8 = readmatrix("../Chaste_connectivity/CC_B3_med_8.csv");
B3_9 = readmatrix("../Chaste_connectivity/CC_B3_med_9.csv");

B3_boundary = [B3_0(:,1);B3_1(:,1);B3_2(:,1);B3_3(:,1);B3_4(:,1);
               B3_5(:,1);B3_6(:,1);B3_7(:,1);B3_8(:,1);B3_9(:,1)];
           
B3_nonboundary = [B3_0(:,3);B3_1(:,3);B3_2(:,3);B3_3(:,3);B3_4(:,3);
               B3_5(:,3);B3_6(:,3);B3_7(:,3);B3_8(:,3);B3_9(:,3) ];
           
B3_curvature = [B3_0(:,2);B3_1(:,2);B3_2(:,2);B3_3(:,2);B3_4(:,2);
               B3_5(:,2);B3_6(:,2);B3_7(:,2);B3_8(:,2);B3_9(:,2) ];



B4_0 = readmatrix("../Chaste_connectivity/CC_B4_med_0.csv");
B4_1 = readmatrix("../Chaste_connectivity/CC_B4_med_1.csv");
B4_2 = readmatrix("../Chaste_connectivity/CC_B4_med_2.csv");
B4_3 = readmatrix("../Chaste_connectivity/CC_B4_med_3.csv");
B4_4 = readmatrix("../Chaste_connectivity/CC_B4_med_4.csv");

B4_5 = readmatrix("../Chaste_connectivity/CC_B4_med_5.csv");
B4_6 = readmatrix("../Chaste_connectivity/CC_B4_med_6.csv");
B4_7 = readmatrix("../Chaste_connectivity/CC_B4_med_7.csv");
B4_8 = readmatrix("../Chaste_connectivity/CC_B4_med_8.csv");
B4_9 = readmatrix("../Chaste_connectivity/CC_B4_med_9.csv");


B4_boundary = [B4_0(:,1);B4_1(:,1);B4_2(:,1);B4_3(:,1);B4_4(:,1);
               B4_5(:,1);B4_6(:,1);B4_7(:,1);B4_8(:,1);B4_9(:,1)];
            
B4_nonboundary = [B4_0(:,3);B4_1(:,3);B4_2(:,3);B4_3(:,3);B4_4(:,3);
               B4_5(:,3);B4_6(:,3);B4_7(:,3);B4_8(:,3);B4_9(:,3) ];
           
B4_curvature = [B4_0(:,2);B4_1(:,2);B4_2(:,2);B4_3(:,2);B4_4(:,2);
               B4_5(:,2);B4_6(:,2);B4_7(:,2);B4_8(:,2);B4_9(:,2) ];

unique(B4_nonboundary)


figure;
subplot(1,4,1)
scatter(B1_boundary./B1_nonboundary,B1_curvature);
ylabel("$\kappa$")
xlabel("$n_{1}/n_{2}$")
title("B1")
ylim([-0.2,0.2])

subplot(1,4,2)
scatter(B2_boundary./B2_nonboundary,B2_curvature);
ylabel("$\kappa$")
xlabel("$n_{1}/n_{2}$")
title("B2")
ylim([-0.2,0.2])

subplot(1,4,3)
scatter(B3_boundary./B3_nonboundary,B3_curvature);
ylabel("$\kappa$")
xlabel("$n_{1}/n_{2}$")
title("B3")
ylim([-0.2,0.2])

subplot(1,4,4)
scatter(B4_boundary./B4_nonboundary,B4_curvature);
ylabel("$\kappa$")
xlabel("$n_{1}/n_{2}$")
title("B4")
ylim([-0.2,0.2])


%
%all_bound = [B1_boundary;B2_boundary;B3_boundary;B4_boundary];
%all_nonbound = [B1_nonboundary;B2_nonboundary;B3_nonboundary;B4_nonboundary];
%all_curv = [B1_curvature;B2_curvature;B3_curvature;B4_curvature];
%%


%B1
n2_1_B1_curv = B1_curvature(B1_nonboundary(:) == 1);
n2_2_B1_curv = B1_curvature(B1_nonboundary(:) == 2);
n2_3_B1_curv = B1_curvature(B1_nonboundary(:) == 3);
n2_4_B1_curv = B1_curvature(B1_nonboundary(:) == 4);


tot_b1 = length(n2_1_B1_curv) + length(n2_2_B1_curv) + length(n2_3_B1_curv) + length(n2_4_B1_curv);
n2_1_B1_pro = length(n2_1_B1_curv)/tot_b1;
n2_2_B1_pro = length(n2_2_B1_curv)/tot_b1;
n2_3_B1_pro = length(n2_3_B1_curv)/tot_b1;
n2_4_B1_pro = length(n2_4_B1_curv)/tot_b1;

n2_1_B1_curv  = [ n2_1_B1_curv; NaN(length(n2_2_B1_curv)-length(n2_1_B1_curv),1)];
n2_2_B1_curv  = [ n2_2_B1_curv; NaN(length(n2_2_B1_curv)-length(n2_2_B1_curv),1)];
n2_3_B1_curv  = [ n2_3_B1_curv; NaN(length(n2_2_B1_curv)-length(n2_3_B1_curv),1)];
n2_4_B1_curv  = [ n2_4_B1_curv; NaN(length(n2_2_B1_curv)-length(n2_4_B1_curv),1)];


%B2
n2_1_B2_curv = B2_curvature(B2_nonboundary(:) == 1);
n2_2_B2_curv = B2_curvature(B2_nonboundary(:) == 2);
n2_3_B2_curv = B2_curvature(B2_nonboundary(:) == 3);
n2_4_B2_curv = B2_curvature(B2_nonboundary(:) == 4);

tot_b2 = length(n2_1_B2_curv) + length(n2_2_B2_curv) + length(n2_3_B2_curv) + length(n2_4_B2_curv);
n2_1_B2_pro = length(n2_1_B2_curv)/tot_b2;
n2_2_B2_pro = length(n2_2_B2_curv)/tot_b2;
n2_3_B2_pro = length(n2_3_B2_curv)/tot_b2;
n2_4_B2_pro = length(n2_4_B2_curv)/tot_b2;

n2_1_B2_curv  = [ n2_1_B2_curv; NaN(length(n2_2_B2_curv)-length(n2_1_B2_curv),1)];
n2_2_B2_curv  = [ n2_2_B2_curv; NaN(length(n2_2_B2_curv)-length(n2_2_B2_curv),1)];
n2_3_B2_curv  = [ n2_3_B2_curv; NaN(length(n2_2_B2_curv)-length(n2_3_B2_curv),1)];
n2_4_B2_curv  = [ n2_4_B2_curv; NaN(length(n2_2_B2_curv)-length(n2_4_B2_curv),1)];

%B3
n2_1_B3_curv = B3_curvature(B3_nonboundary(:) == 1);
n2_2_B3_curv = B3_curvature(B3_nonboundary(:) == 2);
n2_3_B3_curv = B3_curvature(B3_nonboundary(:) == 3);
n2_4_B3_curv = B3_curvature(B3_nonboundary(:) == 4);

tot_b3 = length(n2_1_B3_curv) + length(n2_2_B3_curv) + length(n2_3_B3_curv) + length(n2_4_B3_curv);
n2_1_B3_pro = length(n2_1_B3_curv)/tot_b3;
n2_2_B3_pro = length(n2_2_B3_curv)/tot_b3;
n2_3_B3_pro = length(n2_3_B3_curv)/tot_b3;
n2_4_B3_pro = length(n2_4_B3_curv)/tot_b3;

n2_1_B3_curv  = [ n2_1_B3_curv; NaN(length(n2_2_B3_curv)-length(n2_1_B3_curv),1)];
n2_2_B3_curv  = [ n2_2_B3_curv; NaN(length(n2_2_B3_curv)-length(n2_2_B3_curv),1)];
n2_3_B3_curv  = [ n2_3_B3_curv; NaN(length(n2_2_B3_curv)-length(n2_3_B3_curv),1)];
n2_4_B3_curv  = [ n2_4_B3_curv; NaN(length(n2_2_B3_curv)-length(n2_4_B3_curv),1)];

%B4
n2_1_B4_curv = B4_curvature(B4_nonboundary(:) == 1);
n2_2_B4_curv = B4_curvature(B4_nonboundary(:) == 2);
n2_3_B4_curv = B4_curvature(B4_nonboundary(:) == 3);
n2_4_B4_curv = B4_curvature(B4_nonboundary(:) == 4);


[p,h] = ranksum(n2_1_B4_curv,n2_2_B4_curv)
[p,h] = ranksum(n2_2_B4_curv,n2_3_B4_curv)
[p,h] = ranksum(n2_3_B4_curv,n2_4_B4_curv)


tot_b4 = length(n2_1_B4_curv) + length(n2_2_B4_curv) + length(n2_3_B4_curv) + length(n2_4_B4_curv);
n2_1_B4_pro = length(n2_1_B4_curv)/tot_b4;
n2_2_B4_pro = length(n2_2_B4_curv)/tot_b4;
n2_3_B4_pro = length(n2_3_B4_curv)/tot_b4;
n2_4_B4_pro = length(n2_4_B4_curv)/tot_b4;

n2_1_B4_curv  = [ n2_1_B4_curv; NaN(length(n2_2_B4_curv)-length(n2_1_B4_curv),1)];
n2_2_B4_curv  = [ n2_2_B4_curv; NaN(length(n2_2_B4_curv)-length(n2_2_B4_curv),1)];
n2_3_B4_curv  = [ n2_3_B4_curv; NaN(length(n2_2_B4_curv)-length(n2_3_B4_curv),1)];
n2_4_B4_curv  = [ n2_4_B4_curv; NaN(length(n2_2_B4_curv)-length(n2_4_B4_curv),1)];



figure('Renderer', 'painters', 'Position', [10 10 0.12 0.2])
violinplot([n2_1_B1_curv,n2_2_B1_curv,n2_3_B1_curv,n2_4_B1_curv],["1","2","3","4"]);

ylim([-0.3,0.3])
xlim([0,5])
   box on
   title("B1")
   
   
figure('Renderer', 'painters', 'Position', [10 10 0.12 0.20])
violinplot([n2_1_B2_curv,n2_2_B2_curv,n2_3_B2_curv,n2_4_B2_curv],["1","2","3","4"]);
ylim([-0.3,0.3])
xlim([0,5])
box on
title("B2")

figure('Renderer', 'painters', 'Position', [10 10 0.12 0.2])
violinplot([n2_1_B3_curv,n2_2_B3_curv,n2_3_B3_curv,n2_4_B3_curv],["1","2","3","4"]);
ylim([-0.3,0.3])
xlim([0,5])
box on
title("B3")

figure('Renderer', 'painters', 'Position', [10 10 0.12 0.2])
violinplot([n2_1_B4_curv,n2_2_B4_curv,n2_3_B4_curv,n2_4_B4_curv],["1","2","3","4"]);
ylim([-0.3,0.3])
xlim([0,5])
box on

title("B4")



figure('Renderer', 'painters', 'Position', [10 10 0.12 0.2])
b = bar([n2_1_B1_pro,n2_2_B1_pro ,n2_3_B1_pro ,0],'linewidth',1.5);
b.FaceColor = 'flat';
b.CData(1,:) = [3/255,124/255,162/255];
b.CData(2,:) = [136/255,77/255,134/255];
b.CData(3,:) = [162/255,54/255,76/255];
b.CData(4,:) = [206/255,157/255,69/255];
ylim([0,0.7])
xlim([0,5])
   box on
   title("B1")

   
   figure('Renderer', 'painters', 'Position', [10 10 0.12 0.2])
b = bar([n2_1_B2_pro,n2_2_B2_pro ,n2_3_B2_pro ,n2_4_B2_pro],'linewidth',1.5);
b.FaceColor = 'flat';
b.CData(1,:) = [3/255,124/255,162/255];
b.CData(2,:) = [136/255,77/255,134/255];
b.CData(3,:) = [162/255,54/255,76/255];
b.CData(4,:) = [206/255,157/255,69/255];
ylim([0,0.7])
xlim([0,5])
   box on
   title("B2")
   
   figure('Renderer', 'painters', 'Position', [10 10 0.12 0.2])
b = bar([n2_1_B3_pro,n2_2_B3_pro ,n2_3_B3_pro ,n2_4_B3_pro],'linewidth',1.5);
b.FaceColor = 'flat';
b.CData(1,:) = [3/255,124/255,162/255];
b.CData(2,:) = [136/255,77/255,134/255];
b.CData(3,:) = [162/255,54/255,76/255];
b.CData(4,:) = [206/255,157/255,69/255];
ylim([0,0.7])
xlim([0,5])
   box on
   title("B3")
   
   figure('Renderer', 'painters', 'Position', [10 10 0.12 0.2])
b = bar([n2_1_B4_pro,n2_2_B4_pro ,n2_3_B4_pro ,n2_4_B4_pro],'linewidth',1.5);
b.FaceColor = 'flat';
b.CData(1,:) = [3/255,124/255,162/255];
b.CData(2,:) = [136/255,77/255,134/255];
b.CData(3,:) = [162/255,54/255,76/255];
b.CData(4,:) = [206/255,157/255,69/255];
ylim([0,0.7])
xlim([0,5])
   box on
   title("B4")
   
   

%%
close all
clc
X4 = [n2_1_B4_curv,n2_2_B4_curv,n2_3_B4_curv,n2_4_B4_curv];

[p,h] = ranksum(n2_1_B4_curv,n2_2_B4_curv)
[p,h] = ranksum(n2_2_B4_curv,n2_3_B4_curv)
[p,h] = ranksum(n2_3_B4_curv,n2_4_B4_curv)


%[p,h,stats] = ranksum(n2_1_B4_curv',n2_4_B4_curv')
