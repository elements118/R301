
%Define risk matrix 1-3 (1 low - 2 medium - 3 severe)
C=[1 2 3 3 3;
    1 2 2 3 3;
    1 1 2 2 3;
    1 1 2 2 2;
    1 1 1 1 2];

%Values for each cell
V=[1 2 3 4 5; 6 7 8 9 10; 11 12 13 14 15; 16 17 18 19 20; 21 22 23 24 25];
ax=axes; 
hold on

%Draw colored grid
[X,Y]=meshgrid(1:size(C,1),1:size(C,2));
h=imagesc(X(:),Y(:),flipud(C));

%Define colors (green, yellow, red)
cmap=[0 1 0;1 1 0;1 0 0];

% Add cell values
str = sprintfc('%d',V(:));
text(X(:),Y(:),str);

% Define row and column labels
RowLabels={'row1','row2','row3','row4','row5'};
ColLabels={'col1','col2','col3','col4','col5'};

% Some axes settings
set(gca,'xtick',unique(X),...
    'ytick',unique(Y),...
    'yticklabels',RowLabels,...
    'xticklabels',ColLabels);

% Add colorbar and describe colors
cb=colorbar(ax,'location','southoutside');
set(cb,'ticks',[1.3 2 2.7],...
    'ticklabels',{'Low','Medium','High'},...
    'ticklength',0);
box on
set(ax,'layer','top');
colormap(cmap);
xlabel("Consequence");
ylabel("Likelihood");