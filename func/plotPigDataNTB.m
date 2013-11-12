function plotPigDataNTB(pigDataNTB, hdr)
ncol = size(pigDataNTB,2);
csch = colourscheme;
lw = 2;
f = figure('Units', 'centimeters', 'Position', [0, 0, 29.7, 21], 'Color', 'w');

syscols = [3,4,5];
poxcols = [8,9];
nirshbcols = [10,11,12];

usenirs = ~all(pigDataNTB(:,poxcols) == -1);

cols = length(syscols);
if usenirs 
    cols = cols + length(poxcols);
end
cols = cols + 2;
colours = 'brkg';
a=1;
for n = syscols
    
    subplot(cols,1,a)
    plot(pigDataNTB(:,1), pigDataNTB(:,n), 'color', csch(n).colour, 'LineWidth', lw);
    title(hdr(n))
    %xlabel('elapsed time (s)');
    ymin = min(pigDataNTB(:,n))
    ymax = max(pigDataNTB(:,n));
    if ~(ymin==ymax)
        ylim([min(pigDataNTB(:,n)), max(pigDataNTB(:,n))]);
    end
    set(gca, 'Box', 'off');
    a=a+1;
end
b=1;
if usenirs
    for n = poxcols
        subplot(cols,1,a)
        plot(pigDataNTB(:,1), pigDataNTB(:,n), 'color', csch(n).colour, 'LineWidth', lw);
        title(hdr(n))
        %xlabel('elapsed time (s)');
        ymin = min(pigDataNTB(:,n))
        ymax = max(pigDataNTB(:,n));
        if ~(ymin==ymax)
            ylim([min(pigDataNTB(:,n)), max(pigDataNTB(:,n))]);
        end
        set(gca, 'Box', 'off');
        a=a+1;
        b=b+1;
    end
end

subplot(cols,1,a)
hold on;
title('Hb')
for n = nirshbcols
    plot(pigDataNTB(:,1), pigDataNTB(:,n), 'color', colours(b), 'LineWidth', lw);
    b=b+1;
end
set(gca, 'Box', 'off');

n=13;
a=a+1;
subplot(cols,1,a)
plot(pigDataNTB(:,1), pigDataNTB(:,n), 'color', colours(b), 'LineWidth', lw);
title('CtOx')
set(gca, 'Box', 'off');

spaceplots(f);

    

