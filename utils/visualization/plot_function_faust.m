function plot_function_faust(S1, f)
    X1 = S1.Pts(:,1);
    Y1 = S1.Pts(:,2);
    Z1 = S1.Pts(:,3);
    
    trimesh(S1.surface.TRIV, X1, Y1, Z1, f, 'FaceColor','interp', ...
        'EdgeColor', [0.2 0.2 0.2]); axis equal;
    view([0 90]);
    axis tight
    set(gca,'xtick',[],'ytick',[],'layer','bottom','box','off')
    axis off;
end