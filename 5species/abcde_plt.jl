using DelimitedFiles
using Plots
using StatsPlots
using Colors
using LaTeXStrings


function asym_plots()
	paths = [L"G_{5}",L"G_{5,BC}",L"G_{5,BD}",L"G_{5,CD}",L"G_{5,B}",L"G_{5,C}",L"G_{5,D}"];
	blue = colorant"rgb(38,117,147)";
	orange = colorant"rgb(255,121,23)";

	costs = readdlm("data/asym_costs.dat",',';skipstart=1);

	p = plot(
		yrange=(0,8.5),
		yticks=0:0.5:8.5,
		tickfontsize=10,
		xlabel="pathways",
		xlabelfontsize=10,
		ylabel=L"\chi{(G)}",
		ylabelfontsize=12,
		legendfontsize=12
	);
	p = groupedbar!(
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
	savefig(p,"plots/abcde_cost.pdf");

	return 0
end

