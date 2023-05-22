function fh=venn_plot2(in1,in2)
label_plot='count';
% fh=figure('Color','w','Position',[32,32,255,185]);

%Finding the intersections of the arrays%
Plotting_Interval = 0.01;
Angles_In_Radians = (0: Plotting_Interval: 2*pi);
Circle_Plot = @(X_Offset,Y_Offset,Radius) plot(X_Offset + Radius*cos(Angles_In_Radians),Y_Offset + Radius*sin(Angles_In_Radians));
hold on

%Plotting the 3 circles%
X_Offset_1 = 0; 
Y_Offset_1 = 0; 
Radius_1 = 3;
Circle_1 = Circle_Plot(X_Offset_1,Y_Offset_1,Radius_1);
fill(Circle_1.XData, Circle_1.YData,'r','FaceAlpha',0.2,'LineWidth',1);

X_Offset_2 = 4; 
Y_Offset_2 = 0; 
Radius_2 = 3;
Circle_2 = Circle_Plot(X_Offset_2,Y_Offset_2,Radius_2);
fill(Circle_2.XData, Circle_2.YData,'g','FaceAlpha',0.2,'LineWidth',1);


%Writing all the labels%
if strcmp(label_plot,'item')
    text(X_Offset_1,Y_Offset_1,strjoin(string(in1)),'color','r');
    text(X_Offset_2,Y_Offset_2,strjoin(string(in2)),'color','g');
    text(-1.2,0,strjoin(string(intersect(in1,in2))));
elseif strcmp(label_plot,'proportion')
    text(X_Offset_1,Y_Offset_1,strjoin(string(in1)),'color','r');
    text(X_Offset_2,Y_Offset_2,strjoin(string(in2)),'color','g');
    text(-1.2,0,strjoin(string(intersect(in1,in2))));
elseif strcmp(label_plot,'count')
    text(X_Offset_1,Y_Offset_1,num2str(length(in1)-length(intersect(in1,in2))),'color','r');
    text(X_Offset_2,Y_Offset_2,num2str(length(in2)-length(intersect(in1,in2))),'color','g');
    text((X_Offset_1+X_Offset_2)/2,(Y_Offset_1+Y_Offset_2)/2,num2str(length(intersect(in1,in2))));
end

%Setting the labels to be relative to the centres%
set(findall(gcf,'type','text'),'HorizontalAlignment','center');
axis equal
axis off
end

% 
% A = {'A','B','C','D','E','F','G','H','I','L'};
% B = {'A','B','C','D','E','F','G',};
% C = {'E','F','G','H','I','L'};
% 
% Finding the intersections of the arrays%
% Intersection_AB = intersect(A,B);
% fprintf("Intersection AB: ")
% disp(Intersection_AB);
% fprintf("\n");
% 
% Intersection_BC = intersect(B,C);
% fprintf("Intersection BC: ")
% disp(Intersection_BC);
% fprintf("\n");
% 
% Intersection_AC = intersect(A,C);
% fprintf("Intersection AC: ")
% disp(Intersection_AC);
% fprintf("\n");
% 
% Intersection_ABC = intersect(Intersection_AB,C);
% fprintf("Intersection ABC: ")
% disp(Intersection_ABC);
% fprintf("\n");
% clc;
% 
% clf;
% Plotting_Interval = 0.01;
% Angles_In_Radians = (0: Plotting_Interval: 2*pi);
% Circle_Plot = @(X_Offset,Y_Offset,Radius) plot(X_Offset + Radius*cos(Angles_In_Radians),Y_Offset + Radius*sin(Angles_In_Radians));
% 
% hold on
% Plotting the 3 circles%
% X_Offset_A = 0; Y_Offset_A = 2; Radius_A = 3;
% Circle_A = Circle_Plot(X_Offset_A,Y_Offset_A,Radius_A);
% fill(Circle_A.XData, Circle_A.YData,'r','FaceAlpha',0.2,'LineWidth',1);
% 
% X_Offset_B = -2; Y_Offset_B = -2; Radius_B = 3;
% Circle_B = Circle_Plot(X_Offset_B,Y_Offset_B,Radius_B);
% fill(Circle_B.XData, Circle_B.YData,'g','FaceAlpha',0.2,'LineWidth',1);
% 
% X_Offset_C = 2; Y_Offset_C = -2; Radius_C = 3;
% Circle_Plot(X_Offset_C,Y_Offset_C,Radius_C);
% Circle_C = Circle_Plot(X_Offset_C,Y_Offset_C,Radius_C);
% fill(Circle_C.XData, Circle_C.YData,'b','FaceAlpha',0.2,'LineWidth',1);
% title("Venn Diagram");
% 
% Writing all the labels%
% A_Label = strjoin(string(A));
% text(X_Offset_A,Y_Offset_A,A_Label,'color','r');
% 
% B_Label = strjoin(string(B));
% text(X_Offset_B,Y_Offset_B,B_Label,'color','g');
% 
% C_Label = strjoin(string(C));
% text(X_Offset_C,Y_Offset_C,C_Label,'color','b');
% 
% AB_Label = strjoin(string(Intersection_AB));
% text(-1.2,0,AB_Label);
% 
% BC_Label = strjoin(string(Intersection_BC));
% text(0,-2,BC_Label);
% 
% AC_Label = strjoin(string(Intersection_AC));
% text(1.2,0,AC_Label);
% 
% ABC_Label = strjoin(string(Intersection_ABC));
% text(0,0,ABC_Label);
% 
% Setting the labels to be relative to the centres%
% set(findall(gcf,'type','text'),'HorizontalAlignment','center');
% axis equal
% axis off