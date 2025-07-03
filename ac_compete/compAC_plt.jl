using DelimitedFiles
using Plots
using StatsPlots
using Colors
using LaTeXStrings


function compAC0_plt()
	species = [L"F",L"A",L"B",L"C_1",L"C_2",L"D"];
	paths = [L"G_{AC}",L"G_{AC,1}",L"G_{AC,2}"];
	red = colorant"rgb(226,58,52)";
	orange = colorant"rgb(255,121,23)";
	blue = colorant"rgb(38,117,147)";
	purple = colorant"rgb(113,46,103)";

	costs = readdlm("data/compAC0_costs.dat",',';skipstart=1);

	conc = readdlm("data/compAC0_q.dat",',';skipstart=1);
	conc1 = readdlm("data/compAC1_q.dat",',';skipstart=1);
	diff = log.(conc[2:end,:]'./conc1[1,:]);

	p1 = plot(
		# size=(800,200),
		# left_margin=5Plots.mm,
		# bottom_margin=5Plots.mm,
		yrange=(-3.0,3.0),
		yticks=-3.0:1.0:3.0,
		tickfontsize=10,
		xlabel="species",
		xlabelfontsize=12,
		ylabel=L"\mu_{s}",
		ylabelfontsize=14,
		legendfontsize=12
	);
	hline!([0.0],lc=:black,ls=:dot,label="");
	p1 = groupedbar!(
		species,
		diff,
		color=[purple blue red],
		label=[L"G_{AC}" L"G_{AC,1}" L"G_{AC,2}"],
		foreground_color_legend = nothing,
		background_color_legend = nothing,
		# extra_kwargs=Dict(:subplot=>Dict(:legend_hfactor=>1.5))
	);

	p2 = plot(
		yrange=(0,10),
		yticks=0:1:10,
		tickfontsize=10,
		xlabel="pathways",
		xlabelfontsize=10,
		ylabel=L"\chi(G)",
		ylabelfontsize=12,
		legend=(position=(0.12,0.95)),
		legendfontsize=12
	);
	p2 = groupedbar!(
		paths,
		reverse(costs[:,1:2],dims=2),
		bar_position = :stack,
		color = [orange blue],
		label = [L"\dot{\Delta}\,(G)" L"\dot{\Sigma}\,(G)"],
		foreground_color_legend = nothing,
		background_color_legend = nothing,
		extra_kwargs=Dict(:subplot=>Dict(:legend_hfactor=>1.5))
	);

	if !isdir("plots")
		mkdir("plots");
	end
	savefig(p1,"plots/ACcomp_conc2.pdf");
	savefig(p2,"plots/ACcomp_cost2.pdf");

	return 0
end


function compAC1_plt()
	species = [L"F",L"A",L"B",L"C_1",L"C_2",L"D"];
	paths = [L"G_{AC}",L"G_{AC,1}",L"G_{AC,2}"];
	red = colorant"rgb(226,58,52)";
	orange = colorant"rgb(255,121,23)";
	blue = colorant"rgb(38,117,147)";
	purple = colorant"rgb(113,46,103)";

	costs = readdlm("data/compAC1_costs.dat",',';skipstart=1);

	conc = readdlm("data/compAC1_q.dat",',';skipstart=1);
	diff = log.(conc[2:end,:]'./conc[1,:]);

	p1 = plot(
		# size=(800,200),
		# left_margin=5Plots.mm,
		# bottom_margin=5Plots.mm,
		yrange=(-1.5,1.5),
		yticks=-1.5:0.5:1.5,
		tickfontsize=10,
		xlabel="species",
		xlabelfontsize=12,
		ylabel=L"\mu_{s}",
		ylabelfontsize=14,
		legendfontsize=12
	);
	hline!([0.0],lc=:black,ls=:dot,label="");
	p1 = groupedbar!(
		species,
		diff,
		color=[purple blue red],
		label=[L"G_{AC}" L"G_{AC,1}" L"G_{AC,2}"],
		foreground_color_legend = nothing,
		background_color_legend = nothing,
		# extra_kwargs=Dict(:subplot=>Dict(:legend_hfactor=>1.5))
	);

	p2 = plot(
		yrange=(0,4),
		yticks=0:1:4,
		tickfontsize=10,
		xlabel="pathways",
		xlabelfontsize=10,
		ylabel=L"\chi(G)",
		ylabelfontsize=12,
		legendfontsize=12
	);
	p2 = groupedbar!(
		paths,
		reverse(costs[:,1:2],dims=2),
		bar_position = :stack,
		color = [orange blue],
		label = [L"\dot{\Delta}\,(G)" L"\dot{\Sigma}\,(G)"],
		foreground_color_legend = nothing,
		background_color_legend = nothing,
		extra_kwargs=Dict(:subplot=>Dict(:legend_hfactor=>1.5))
	);

	if !isdir("plots")
		mkdir("plots");
	end
	savefig(p1,"plots/ACcomp_conc1.pdf");
	savefig(p2,"plots/ACcomp_cost1.pdf");

	return 0
end

