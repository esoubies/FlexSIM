function PlotSimuRes(params)
%--------------------------------------------------------------------------
% Function params = PlotResults(params)
% 
% Loads the table stored at `Results.m` with accumulated results of
% parameter estimation through different methodologies. Plots the absolute
% error in the wavevector estmation, in the phase estimation and in the
% Frobenius norm of the reconstructed patterns. Does not need any innputs
% nor generate any outputs.
%
%
% Copyright (2022) A. Nogueron (anogueron.1996@gmail.com)
%                  E. Soubies (emmanuel.soubies@irit.fr) 
%--------------------------------------------------------------------------
T = load('Results.mat').T; 
% Generate errors for wavevectors
T{:, "K PattErr"} = vecnorm(T{:, "K (GT)"} - T{:, "K (From Patt)"}, 2, 2);
T{:, "K CC eq-ph Err"} = vecnorm(T{:, "K (GT)"} - T{:, "K (CC eq-ph)"}, 2, 2);
T{:, "K CC eq-ph Ref Err"} = vecnorm(T{:, "K (GT)"} - T{:, "K (CC eq-ph Ref)"}, 2, 2);
T{:, "K CC eq-ph Filt Err"} = vecnorm(T{:, "K (GT)"} - T{:, "K (CC eq-ph Filt)"}, 2, 2);
T{:, "K CC Err"} = vecnorm(T{:, "K (GT)"} - T{:, "K (CC)"}, 2, 2);
T{:, "K CC Ref Err"} = vecnorm(T{:, "K (GT)"} - T{:, "K (CC Ref)"}, 2, 2);
T{:, "K CC Filt Err"} = vecnorm(T{:, "K (GT)"} - T{:, "K (CC Filt)"}, 2, 2);

T{:, "Ph #1 PattErr"} = vecnorm(T{:, "Ph #1 (GT)"} - T{:, "Ph #1 (from Patt)"}, 2, 2);
T{:, "Ph #1 CC eq-ph Err"} = vecnorm(T{:, "Ph #1 (GT)"} - T{:, "Ph #1 (CC eq-ph)"}, 2, 2);
T{:, "Ph #1 CC eq-ph Ref Err"} = vecnorm(T{:, "Ph #1 (GT)"} - T{:, "Ph #1 (CC eq-ph Ref)"}, 2, 2);
T{:, "Ph #1 CC eq-ph Filt Err"} = vecnorm(T{:, "Ph #1 (GT)"} - T{:, "Ph #1 (CC eq-ph Filt)"}, 2, 2);
T{:, "Ph #1 CC Err"} = vecnorm(T{:, "Ph #1 (GT)"} - T{:, "Ph #1 (CC)"}, 2, 2);
T{:, "Ph #1 CC Ref Err"} = vecnorm(T{:, "Ph #1 (GT)"} - T{:, "Ph #1 (CC Ref)"}, 2, 2);
T{:, "Ph #1 CC Filt Err"} = vecnorm(T{:, "Ph #1 (GT)"} - T{:, "Ph #1 (CC Filt)"}, 2, 2);

T{:, "Ph #2 PattErr"} = vecnorm(T{:, "Ph #2 (GT)"} - T{:, "Ph #2 (from Patt)"}, 2, 2);
T{:, "Ph #2 CC eq-ph Err"} = vecnorm(T{:, "Ph #2 (GT)"} - T{:, "Ph #2 (CC eq-ph)"}, 2, 2);
T{:, "Ph #2 CC eq-ph Ref Err"} = vecnorm(T{:, "Ph #2 (GT)"} - T{:, "Ph #2 (CC eq-ph Ref)"}, 2, 2);
T{:, "Ph #2 CC eq-ph Filt Err"} = vecnorm(T{:, "Ph #2 (GT)"} - T{:, "Ph #2 (CC eq-ph Filt)"}, 2, 2);
T{:, "Ph #2 CC Err"} = vecnorm(T{:, "Ph #2 (GT)"} - T{:, "Ph #2 (CC)"}, 2, 2);
T{:, "Ph #2 CC Ref Err"} = vecnorm(T{:, "Ph #2 (GT)"} - T{:, "Ph #2 (CC Ref)"}, 2, 2);
T{:, "Ph #2 CC Filt Err"} = vecnorm(T{:, "Ph #2 (GT)"} - T{:, "Ph #2 (CC Filt)"}, 2, 2);

T{:, "Ph #3 PattErr"} = vecnorm(T{:, "Ph #3 (GT)"} - T{:, "Ph #3 (from Patt)"}, 2, 2);
T{:, "Ph #3 CC eq-ph Err"} = vecnorm(T{:, "Ph #3 (GT)"} - T{:, "Ph #3 (CC eq-ph)"}, 2, 2);
T{:, "Ph #3 CC eq-ph Ref Err"} = vecnorm(T{:, "Ph #3 (GT)"} - T{:, "Ph #3 (CC eq-ph Ref)"}, 2, 2);
T{:, "Ph #3 CC eq-ph Filt Err"} = vecnorm(T{:, "Ph #3 (GT)"} - T{:, "Ph #3 (CC eq-ph Filt)"}, 2, 2);
T{:, "Ph #3 CC Err"} = vecnorm(T{:, "Ph #3 (GT)"} - T{:, "Ph #3 (CC)"}, 2, 2);
T{:, "Ph #3 CC Ref Err"} = vecnorm(T{:, "Ph #3 (GT)"} - T{:, "Ph #3 (CC Ref)"}, 2, 2);
T{:, "Ph #3 CC Filt Err"} = vecnorm(T{:, "Ph #3 (GT)"} - T{:, "Ph #3 (CC Filt)"}, 2, 2);

% Calculate initialization pixel resolution in Fourier space (assuming padding factor was 4)
padSz = params.sz(1:2)*(params.fac+1); 
pixRes = pi./params.res./padSz; 

% Extract the existing groups in the table
grps = grpstats(T{:,"K PattErr"}, {T{:,"Contrast"}, T{:,"MEP"}}, "gname"); nbGrps = size(grps, 1); 
figure; sgtitle("Wavevector estimation error"); spCount = 1; 
pause('on') % For plotting purposes
for i = 1:nbGrps
    i
    subplot(3,3,i)
    % Extract group contrast and noise level + declare subtable 
    a = str2double(cell2mat(grps(i, 1))); MEP = str2double(cell2mat(grps(i, 2))); 
    idxGrp = T{:,"Contrast"} == a & T{:,"MEP"} == MEP; subT = T(idxGrp, :);
    
%     figure; subplot(1, 2, 1); 
    % K boxplot
    x = [subT{:, "K PattErr"}; subT{:, "K CC eq-ph Err"}; subT{:, "K CC eq-ph Ref Err"}; subT{:, "K CC eq-ph Filt Err"}; 
        subT{:, "K CC Err"}; subT{:, "K CC Ref Err"}; subT{:, "K CC Filt Err"}];
    g = [repmat({"K PattErr"}, length(subT{:, "K PattErr"}), 1); 
        repmat({"K CC eq-ph Err"}, length(subT{:, "K CC eq-ph Err"}), 1); repmat({"K CC eq-ph Ref Err"}, length(subT{:, "K CC eq-ph Ref Err"}), 1); repmat({"K CC eq-ph Filt Err"}, length(subT{:, "K CC eq-ph Filt Err"}), 1);...
        repmat({"K CC Err"}, length(subT{:, "K CC Err"}), 1); repmat({"K CC Ref Err"}, length(subT{:, "K CC Ref Err"}), 1); repmat({"K CC Filt Err"}, length(subT{:, "K CC Filt Err"}), 1);];
    boxplot(x, g); 
    title(sprintf("MEP: %d, Contrast: %g",  MEP, a), 'Interpreter','latex', 'FontSize', 16, 'FontName','TimesNewRoman');
    if i < 7; set(gca,'XTickLabel',[]); end
    if ismember(i, [1, 4, 7]); ylabel('$\|\mathbf{k}_0 - \mathbf{k}_{est}\|_2$', 'Interpreter','latex', 'FontSize', 16, 'FontName','TimesNewRoman'); end
    grid on    
    h(1) = line([0 7.5], [pi./params.res./padSz(1)/2 pi./params.res./padSz(1)/2], "LineStyle", "--", "Color", "g"); 
    if padSz(1) ~= padSz(2)
        h(2) = line([0 7.5], [pi./params.res./padSz(2)/2 pi./params.res./padSz(2)/2], "LineStyle", "--", "Color", "y"); 
        legend([h(1), h(2)], "Half pixel resolution (x)", "Half pixel resolution (y)")
    else
        legend([h(1)], "Half pixel resolution")
    end
    drawnow; pause(0.1);
end

% Phase (all three) boxplot
figure; sgtitle("Phase estimation error");
for i = 1:nbGrps
    subplot(3,3,i)
    % Extract group contrast and noise level + declare subtable 
    a = str2double(cell2mat(grps(i, 1))); MEP = str2double(cell2mat(grps(i, 2))); 
    idxGrp = T{:,"Contrast"} == a & T{:,"MEP"} == MEP; subT = T(idxGrp, :);

    x = [subT{:, "Ph #1 PattErr"}; subT{:, "Ph #2 PattErr"}; subT{:, "Ph #3 PattErr"};
        subT{:, "Ph #1 CC eq-ph Err"}; subT{:, "Ph #2 CC eq-ph Err"}; subT{:, "Ph #3 CC eq-ph Err"};
        subT{:, "Ph #1 CC eq-ph Ref Err"}; subT{:, "Ph #2 CC eq-ph Ref Err"}; subT{:, "Ph #3 CC eq-ph Ref Err"};
        subT{:, "Ph #1 CC eq-ph Filt Err"}; subT{:, "Ph #2 CC eq-ph Filt Err"}; subT{:, "Ph #3 CC eq-ph Filt Err"};
        subT{:, "Ph #1 CC Err"}; subT{:, "Ph #2 CC Err"}; subT{:, "Ph #3 CC Err"}; 
        subT{:, "Ph #1 CC Ref Err"}; subT{:, "Ph #2 CC Ref Err"}; subT{:, "Ph #3 CC Ref Err"}; 
        subT{:, "Ph #1 CC Filt Err"}; subT{:, "Ph #2 CC Filt Err"}; subT{:, "Ph #3 CC Filt Err"}];
    g = [repmat({"Ph PattErr"}, 3*length(subT{:, "Ph #1 PattErr"}), 1); 
        repmat({"Ph CC eq-ph Err"}, 3*length(subT{:, "Ph #1 CC eq-ph Err"}), 1); repmat({"Ph CC eq-ph Ref Err"}, 3*length(subT{:, "Ph #1 CC eq-ph Ref Err"}), 1); repmat({"Ph CC eq-ph Filt Err"}, 3*length(subT{:, "Ph #1 CC eq-ph Filt Err"}), 1);...
        repmat({"Ph CC Err"}, 3*length(subT{:, "Ph #1 CC Err"}), 1); repmat({"Ph CC Ref Err"}, 3*length(subT{:, "Ph #1 CC Ref Err"}), 1); repmat({"Ph CC Filt Err"}, 3*length(subT{:, "Ph #1 CC Filt Err"}), 1);];
    boxplot(x, g); 
    title(sprintf("MEP: %d, Contrast: %g",  MEP, a), 'Interpreter','latex', 'FontSize', 16, 'FontName','TimesNewRoman');
    if i < 7; set(gca,'XTickLabel',[]); end
    if ismember(i, [1, 4, 7]); ylabel('$\|\phi_0 - \phi_{est}\| [^\circ]$', 'Interpreter','latex', 'FontSize', 16, 'FontName','TimesNewRoman'); end
    grid on; drawnow; pause(1);    
%     sgtitle(sprintf("Parameter Estmation Error (MEP: %d, Contrast: %g)",  MEP, a), 'FontSize', 16, 'FontName','TimesNewRoman');
end
pause('off') % For plotting purposes
% Extract Children to plot median graphs
% hold on
% a = get(get(gca,'children'),'children'); t = get(a,'tag'); idx=strcmpi(t,'median');
% medLine=a(idx); 
% % Extract stats for k detected from pattern
% [mean_kPatt, std_kPatt, grps_kPatt] = grpstats(T{:,"K PattErr"}, {T{:,"Contrast"}, T{:,"MEP"}}, ["mean", "std", "gname"]);
% 
% C = {'k','m','g','y'};
% h = zeros(1, length(medLine)); 
% leg = zeros(1, length(medLine));
% for i = [1:length(mean_kPatt)] % Plot lines for the specific 
%     line(medLine(end).XData, [mean_kPatt(i)]);
%     leg(i) = strcat("Contrast: ",grps_kPatt(i, 1), ", MEP: ", grps_kPatt(i, 2));
% end

% Extracting stats
% [mean_kPatt, std_kPatt, grps_kPatt] = grpstats(T{:,"K PattErr"}, {T{:,"Contrast"}, T{:,"MEP"}}, ["mean", "std", "gname"]);
% % Some dirty coding to get doubles out of the groups info
% grps_kPatt = cell2table(grps_kPatt); grps_kPatt.grps_kPatt1 = str2num(cell2mat(grps_kPatt.grps_kPatt1));
% grps_kPatt.grps_kPatt2 = cell2mat(grps_kPatt.grps_kPatt2); grps_kPatt = table2array(grps_kPatt);
% [mean_kCCEqPh, std_kCCEqPh, grps_kCCEqPh] = grpstats(T{:,"K CC eq-ph Err"}, {T{:,"Contrast"}, T{:,"MEP"}}, ["mean", "std", "gname"]);
% [mean_kCC, std_kCC, grps_kCC] = grpstats(T{:,"K CC Err"}, {T{:,"Contrast"}, T{:,"MEP"}}, ["mean", "std", "gname"]);
% [mean_phPatt, std_phPatt, grps_phPatt] = grpstats(T{:,"Ph #1 PattErr"}, {T{:,"Contrast"}, T{:,"MEP"}}, ["mean", "std", "gname"]);
% [mean_phCCEqPh, std_phCCEqPh, grps_phCCEqPh] = grpstats(T{:,"Ph #1 CC eq-ph Err"}, {T{:,"Contrast"}, T{:,"MEP"}}, ["mean", "std", "gname"]);
% [mean_phCC, std_phCC, grps_phCC] = grpstats(T{:,"Ph #1 CC Err"}, {T{:,"Contrast"}, T{:,"MEP"}}, ["mean", "std", "gname"]);
% [mean_froPatt, std_froPatt, grps_froPatt] = grpstats(T{:,"K PattErr"}, {T{:,"Contrast"}, T{:,"MEP"}}, ["mean", "std", "gname"]);
% [mean_froCCEqPh, std_froCCEqPh, grps_froCCEqPh] = grpstats(T{:,"K CC eq-ph Err"}, {T{:,"Contrast"}, T{:,"MEP"}}, ["mean", "std", "gname"]);
% [mean_froCC, std_froCC, grps_froCC] = grpstats(T{:,"K CC Err"}, {T{:,"Contrast"}, T{:,"MEP"}}, ["mean", "std", "gname"]);


% Actual plotting
% [X, Y] = meshgrid(grps_kPatt(:,1), grps_kPatt(:,2));
% surf(X, Y, mean_kPatt);
end