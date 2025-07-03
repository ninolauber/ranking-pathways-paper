using DelimitedFiles
using Plots
using StatsPlots
using Colors
using LaTeXStrings


function sym_plots()
	species = [L"A",L"B",L"C",L"D"];
	paths = [L"G_{4}",L"G_{4,B}",L"G_{4,C}"];
	red = colorant"rgb(226,58,52)";
	orange = colorant"rgb(255,121,23)";
	blue = colorant"rgb(38,117,147)";
	purple = colorant"rgb(113,46,103)";

	costs = readdlm("data/sym_costs.dat",',';skipstart=1);

	conc = readdlm("data/sym_q.dat",',';skipstart=1);
	diff = log.(conc[2:end,:]'./conc[1,:]);

	p1 = plot(
		size=(800,200),
		left_margin=5Plots.mm,
		bottom_margin=5Plots.mm,
		yrange=(-0.5,2.8),
		yticks=-0.5:0.5:2.8,
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
		label=[L"G_{4}" L"G_{4,B}" L"G_{4,C}"],
		foreground_color_legend = nothing,
		background_color_legend = nothing,
		# extra_kwargs=Dict(:subplot=>Dict(:legend_hfactor=>1.5))
	);

	p2 = plot(
		yrange=(0,4),
		yticks=0:0.5:4,
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
	savefig(p1,"plots/abcd_conc1.pdf");
	savefig(p2,"plots/abcd_cost1.pdf");

	return 0
end


function asym1_plots()
	species = [L"A",L"B",L"C",L"D"];
	paths = [L"G_{4}",L"G_{4,B}",L"G_{4,C}"];
	red = colorant"rgb(226,58,52)";
	orange = colorant"rgb(255,121,23)";
	blue = colorant"rgb(38,117,147)";
	purple = colorant"rgb(113,46,103)";

	costs1 = readdlm("data/asym11_costs.dat",',';skipstart=1);
	costs2 = readdlm("data/asym12_costs.dat",',';skipstart=1);

	conc1 = readdlm("data/asym11_q.dat",',';skipstart=1);
	conc2 = readdlm("data/asym12_q.dat",',';skipstart=1);
	diff1 = log.(conc1[2:end,:]'./conc1[1,:]);
	diff2 = log.(conc2[2:end,:]'./conc2[1,:]);
	diff = [diff1[:,1] diff2[:,1] diff1[:,2] diff2[:,2] diff1[:,3] diff2[:,3]];

	p1 = plot(
		size=(800,200),
		left_margin=5Plots.mm,
		bottom_margin=5Plots.mm,

		yrange=(-0.6,2.8),
		yticks=-1:0.5:3,
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
		color=[purple purple blue blue red red],
		alpha = [1.0 0.5 1.0 0.5 1.0 0.5],		
		label = [L"G_{4}" "" L"G_{4,B}" "" L"G_{4,C}" ""],
		foreground_color_legend = nothing,
		background_color_legend = nothing,
		# extra_kwargs=Dict(:subplot=>Dict(:legend_hfactor=>1.5))
	);

	p2 = plot(
		xrange=(0.5,3.5),
		xticks=(1:1:3,paths),
		yrange=(0,7),
		yticks=0:0.5:7,
		tickfontsize=10,
		xlabel="pathways",
		xlabelfontsize=10,
		ylabel=L"\chi(G)",
		ylabelfontsize=12,
		legendfontsize=12
	);
	p2 = groupedbar!(
		[0.75,1.75,2.75],
		reverse(costs1[:,1:2],dims=2),
		bar_position = :stack,
		bar_width=0.3,
		color = [orange blue],
		alpha = 1.0,
		label = [L"\dot{\Delta}\,(G)" L"\dot{\Sigma}\,(G)"],
		foreground_color_legend = nothing,
		background_color_legend = nothing,
		extra_kwargs=Dict(:subplot=>Dict(:legend_hfactor=>1.5))
	);
	p2 = groupedbar!(
		[1.25,2.25,3.25],
		reverse(costs2[:,1:2],dims=2),
		bar_position = :stack,
		bar_width=0.3,
		color = [orange blue],
		alpha = 0.5,
		label = ["" ""]
	);

	if !isdir("plots")
		mkdir("plots");
	end
	savefig(p1,"plots/abcd_conc2.pdf");
	savefig(p2,"plots/abcd_cost2.pdf");

	return 0
end


function asym2_plots()
	species = [L"A",L"B",L"C",L"D"];
	paths = [L"G_{4}",L"G_{4,B}",L"G_{4,C}"];
	red = colorant"rgb(226,58,52)";
	orange = colorant"rgb(255,121,23)";
	blue = colorant"rgb(38,117,147)";
	purple = colorant"rgb(113,46,103)";

	costs1 = readdlm("data/asym21_costs.dat",',';skipstart=1);
	costs2 = readdlm("data/asym22_costs.dat",',';skipstart=1);

	conc1 = readdlm("data/asym21_q.dat",',';skipstart=1);
	conc2 = readdlm("data/asym22_q.dat",',';skipstart=1);
	diff1 = log.(conc1[2:end,:]'./conc1[1,:]);
	diff2 = log.(conc2[2:end,:]'./conc2[1,:]);
	diff = [diff1[:,1] diff2[:,1] diff1[:,2] diff2[:,2] diff1[:,3] diff2[:,3]];

	p1 = plot(
		size=(800,200),
		left_margin=5Plots.mm,
		bottom_margin=5Plots.mm,

		yrange=(-3.5,4),
		yticks=-4:1:4,
		tickfontsize=10,
		xlabel="species",
		xlabelfontsize=12,
		ylabel=L"\mu_{S}",
		ylabelfontsize=14,
		legendfontsize=8
	);
	hline!([0.0],lc=:black,ls=:dot,label="");
	p1 = groupedbar!(
		species,
		diff,
		color=[purple purple blue blue red red],
		alpha = [1.0 0.5 1.0 0.5 1.0 0.5],
		label = [L"G_{4}" "" L"G_{4,B}" "" L"G_{4,C}" ""],
		foreground_color_legend = nothing,
		background_color_legend = nothing,
		# extra_kwargs=Dict(:subplot=>Dict(:legend_hfactor=>1.5))
	);

	p2 = plot(
		xrange=(0.5,3.5),
		xticks=(1:1:3,paths),
		yrange=(0,13),
		yticks=0:1:13,
		tickfontsize=10,
		xlabel="pathways",
		xlabelfontsize=10,
		ylabel=L"\chi(G)",
		ylabelfontsize=12,
		legendfontsize=12
	);
	p2 = groupedbar!(
		[0.75,1.75,2.75],
		reverse(costs1[:,1:2],dims=2),
		bar_position = :stack,
		bar_width=0.3,
		color = [orange blue],
		alpha = 1.0,
		label = [L"\dot{\Delta}\,(G)" L"\dot{\Sigma}\,(G)"],
		foreground_color_legend = nothing,
		background_color_legend = nothing,
		extra_kwargs=Dict(:subplot=>Dict(:legend_hfactor=>1.5))
	);
	p2 = groupedbar!(
		[1.25,2.25,3.25],
		reverse(costs2[:,1:2],dims=2),
		bar_position = :stack,
		bar_width=0.3,
		color = [orange blue],
		alpha = 0.5,
		label = ["" ""]
	);

	if !isdir("plots")
		mkdir("plots");
	end
	savefig(p1,"plots/abcd_conc3.pdf");
	savefig(p2,"plots/abcd_cost3.pdf");

	return 0
end

