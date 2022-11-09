function PlotSimuRes()
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
% For Phase (To deprecate?)
% T{:, "Ph #1 PattErr"} = vecnorm(T{:, "Ph #1 (GT)"} - T{:, "Ph #1 (from Patt)"}, 2, 2);
% T{:, "Ph #1 CC eq-ph Err"} = vecnorm(T{:, "Ph #1 (GT)"} - T{:, "Ph #1 (CC eq-ph)"}, 2, 2);
% T{:, "Ph #1 CC Err"} = vecnorm(T{:, "Ph #1 (CC)"} - T{:, "Ph #1 (CC)"}, 2, 2);
% T{:, "Ph #1 Ref Err"} = vecnorm(T{:, "Ph #1 (Ref)"} - T{:, "Ph #1 (Ref)"}, 2, 2);
% T{:, "Ph #1 Filt Err"} = vecnorm(T{:, "Ph #1 (Filt)"} - T{:, "Ph #1 (Filt)"}, 2, 2);
% T{:, "Ph #2 PattErr"} = vecnorm(T{:, "Ph #2 (GT)"} - T{:, "Ph #2 (from Patt)"}, 2, 2);
% T{:, "Ph #2 CC eq-ph Err"} = vecnorm(T{:, "Ph #2 (GT)"} - T{:, "Ph #2 (CC eq-ph)"}, 2, 2);
% T{:, "Ph #2 Ref Err"} = vecnorm(T{:, "Ph #2 (Ref)"} - T{:, "Ph #2 (Ref)"}, 2, 2);
% T{:, "Ph #2 Filt Err"} = vecnorm(T{:, "Ph #2 (Filt)"} - T{:, "Ph #2 (Filt)"}, 2, 2);
% T{:, "Ph #2 CC Err"} = vecnorm(T{:, "Ph #2 (CC)"} - T{:, "Ph #2 (CC)"}, 2, 2);
% T{:, "Ph #3 PattErr"} = vecnorm(T{:, "Ph #3 (GT)"} - T{:, "Ph #3 (from Patt)"}, 2, 2);
% T{:, "Ph #3 CC eq-ph Err"} = vecnorm(T{:, "Ph #3 (GT)"} - T{:, "Ph #3 (CC eq-ph)"}, 2, 2);
% T{:, "Ph #3 CC Err"} = vecnorm(T{:, "Ph #3 (CC)"} - T{:, "Ph #3 (CC)"}, 2, 2);
% T{:, "Ph #3 Ref Err"} = vecnorm(T{:, "Ph #3 (Ref)"} - T{:, "Ph #3 (Ref)"}, 2, 2);
% T{:, "Ph #3 Filt Err"} = vecnorm(T{:, "Ph #3 (Filt)"} - T{:, "Ph #3 (Filt)"}, 2, 2);
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

% T{:, "Ph #1 Arg Err"} = vecnorm(T{:, "Ph #1 (GT)"} - T{:, "Ph #1 (Arg)"}, 2, 2);
% T{:, "Ph #1 J Err"} = vecnorm(T{:, "Ph #1 (GT)"} - T{:, "Ph #1 (J)"}, 2, 2);
% T{:, "Ph #1 Filt Err"} = vecnorm(T{:, "Ph #1 (GT)"} - T{:, "Ph #1 (Filt)"}, 2, 2);
% % T{:, "Ph #2 Arg Err"} = vecnorm(T{:, "Ph #2 (GT)"} - T{:, "Ph #2 (Arg)"}, 2, 2);
% T{:, "Ph #2 Filt Err"} = vecnorm(T{:, "Ph #2 (GT)"} - T{:, "Ph #2 (Filt)"}, 2, 2);
% T{:, "Ph #2 J Err"} = vecnorm(T{:, "Ph #2 (GT)"} - T{:, "Ph #2 (J)"}, 2, 2);
% % T{:, "Ph #3 Arg Err"} = vecnorm(T{:, "Ph #3 (GT)"} - T{:, "Ph #3 (Arg)"}, 2, 2);
% T{:, "Ph #3 J Err"} = vecnorm(T{:, "Ph #3 (GT)"} - T{:, "Ph #3 (J)"}, 2, 2);
% T{:, "Ph #3 Filt Err"} = vecnorm(T{:, "Ph #3 (GT)"} - T{:, "Ph #3 (Filt)"}, 2, 2);


% K boxplot
x = [T{:, "K PattErr"}; T{:, "K CC eq-ph Err"}; T{:, "K CC eq-ph Ref Err"}; T{:, "K CC eq-ph Filt Err"}; 
    T{:, "K CC Err"}; T{:, "K CC Ref Err"}; T{:, "K CC Filt Err"}];
g = [repmat({"K PattErr"}, length(T{:, "K PattErr"}), 1); 
    repmat({"K CC eq-ph Err"}, length(T{:, "K CC eq-ph Err"}), 1); repmat({"K CC eq-ph Ref Err"}, length(T{:, "K CC eq-ph Ref Err"}), 1); repmat({"K CC eq-ph Filt Err"}, length(T{:, "K CC eq-ph Filt Err"}), 1);...
    repmat({"K CC Err"}, length(T{:, "K CC Err"}), 1); repmat({"K CC Ref Err"}, length(T{:, "K CC Ref Err"}), 1); repmat({"K CC Filt Err"}, length(T{:, "K CC Filt Err"}), 1);];
figure; boxplot(x, g); title("Boxplot of Wavevector Error", 'Interpreter','latex', 'FontSize', 12, 'FontName','TimesNewRoman');
ylabel('$\|\mathbf{k}_0 - \mathbf{k}_{est}\|_2$', 'Interpreter','latex', 'FontSize', 16, 'FontName','TimesNewRoman');
grid on
ylim([0 10e-5])
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

% Ph #1 boxplot (Obsolete w previous phase methods)
% x = [T{:, "Ph #1 PattErr"}; T{:, "Ph #1 CC eq-ph Err"}; T{:, "Ph #1 CC Err"}; T{:, "Ph #1 Ref Err"}; T{:, "Ph #1 Filt Err"}];
% g = [repmat({"Ph #1 PattErr"}, length(T{:, "Ph #1 PattErr"}), 1); repmat({"Ph #1 CC eq-ph Err"}, length(T{:, "Ph #1 CC eq-ph Err"}), 1);...
%     repmat({"Ph #1 CC Err"}, length(T{:, "Ph #1 CC Err"}), 1); repmat({"Ph #1 Ref Err"}, length(T{:, "Ph #1 Ref Err"}), 1); repmat({"Ph #1 Filt Err"}, length(T{:, "Ph #1 Filt Err"}), 1);]
% boxplot(x, g)

% Ph #1 boxplot (w proposed methods w and wo filter)
% x = rad2deg([T{:, "Ph #1 J Err"}; T{:, "Ph #1 Filt Err"}]);
% g = [repmat({"Ph #1 J Err"}, length(T{:, "Ph #1 J Err"}), 1); repmat({"Ph #1 Filt Err"}, length(T{:, "Ph #1 Filt Err"}), 1);];
x = [T{:, "Ph #1 PattErr"}; T{:, "Ph #1 CC eq-ph Err"}; T{:, "Ph #1 CC eq-ph Ref Err"}; T{:, "Ph #1 CC eq-ph Filt Err"}; 
    T{:, "Ph #1 CC Err"}; T{:, "Ph #1 CC Ref Err"}; T{:, "Ph #1 CC Filt Err"}];
g = [repmat({"Ph #1 PattErr"}, length(T{:, "Ph #1 PattErr"}), 1); 
    repmat({"Ph #1 CC eq-ph Err"}, length(T{:, "Ph #1 CC eq-ph Err"}), 1); repmat({"Ph #1 CC eq-ph Ref Err"}, length(T{:, "Ph #1 CC eq-ph Ref Err"}), 1); repmat({"Ph #1 CC eq-ph Filt Err"}, length(T{:, "Ph #1 CC eq-ph Filt Err"}), 1);...
    repmat({"Ph #1 CC Err"}, length(T{:, "Ph #1 CC Err"}), 1); repmat({"Ph #1 CC Ref Err"}, length(T{:, "Ph #1 CC Ref Err"}), 1); repmat({"Ph #1 CC Filt Err"}, length(T{:, "Ph #1 CC Filt Err"}), 1);];
figure; boxplot(x, g); title("Boxplot of Phase Estimation Error", 'Interpreter','latex', 'FontSize', 12, 'FontName','TimesNewRoman');
ylabel('$\|\phi_0 - \phi_{est}\| [^\circ]$', 'Interpreter','latex', 'FontSize', 16, 'FontName','TimesNewRoman');
grid on

% Extracting stats
[mean_kPatt, std_kPatt, grps_kPatt] = grpstats(T{:,"K PattErr"}, {T{:,"Contrast"}, T{:,"MEP"}}, ["mean", "std", "gname"]);
% Some dirty coding to get doubles out of the groups info
grps_kPatt = cell2table(grps_kPatt); grps_kPatt.grps_kPatt1 = str2num(cell2mat(grps_kPatt.grps_kPatt1));
grps_kPatt.grps_kPatt2 = cell2mat(grps_kPatt.grps_kPatt2); grps_kPatt = table2array(grps_kPatt);
[mean_kCCEqPh, std_kCCEqPh, grps_kCCEqPh] = grpstats(T{:,"K CC eq-ph Err"}, {T{:,"Contrast"}, T{:,"MEP"}}, ["mean", "std", "gname"]);
[mean_kCC, std_kCC, grps_kCC] = grpstats(T{:,"K CC Err"}, {T{:,"Contrast"}, T{:,"MEP"}}, ["mean", "std", "gname"]);
[mean_phPatt, std_phPatt, grps_phPatt] = grpstats(T{:,"Ph #1 PattErr"}, {T{:,"Contrast"}, T{:,"MEP"}}, ["mean", "std", "gname"]);
[mean_phCCEqPh, std_phCCEqPh, grps_phCCEqPh] = grpstats(T{:,"Ph #1 CC eq-ph Err"}, {T{:,"Contrast"}, T{:,"MEP"}}, ["mean", "std", "gname"]);
[mean_phCC, std_phCC, grps_phCC] = grpstats(T{:,"Ph #1 CC Err"}, {T{:,"Contrast"}, T{:,"MEP"}}, ["mean", "std", "gname"]);
[mean_froPatt, std_froPatt, grps_froPatt] = grpstats(T{:,"K PattErr"}, {T{:,"Contrast"}, T{:,"MEP"}}, ["mean", "std", "gname"]);
[mean_froCCEqPh, std_froCCEqPh, grps_froCCEqPh] = grpstats(T{:,"K CC eq-ph Err"}, {T{:,"Contrast"}, T{:,"MEP"}}, ["mean", "std", "gname"]);
[mean_froCC, std_froCC, grps_froCC] = grpstats(T{:,"K CC Err"}, {T{:,"Contrast"}, T{:,"MEP"}}, ["mean", "std", "gname"]);


% Actual plotting
% [X, Y] = meshgrid(grps_kPatt(:,1), grps_kPatt(:,2));
% surf(X, Y, mean_kPatt);
end