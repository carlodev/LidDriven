using JLD2
using DataFrames, XLSX, CSV
using Plots
using Statistics

#Re 1000
Bmk_1000 = DataFrame(XLSX.readtable("Benchmark/LidDriven_ux_B.xlsx", "Re1000")...)
FEM_1000 = DataFrame(XLSX.readtable("Results/Re1000_ux.xlsx", "Sheet1")...)
plot(title="Lid Driven Cavity Flow, Re=1000", Bmk_1000.y, Bmk_1000.x0p5, seriestype = :scatter, label="Benchmark")
xlabel!("y, m")
ylabel!("Velocity in x direction, m/s")
plot!(FEM_1000.arc_length,FEM_1000.uhx, linestyle=:dot, label="FEM")
plot!(legend=:bottomright)
savefig("LidDriven_Re_1000.pdf")


#Re 5000
Bmk_5000 = DataFrame(XLSX.readtable("Benchmark/LidDriven_ux_B.xlsx", "Re5000")...)
FEM_5000 = DataFrame(XLSX.readtable("Results/Re5000_ux.xlsx", "Sheet1")...)
plot(title="Lid Driven Cavity Flow, Re=5000", Bmk_5000.y, Bmk_5000.x0p5, seriestype = :scatter, label="Benchmark")
xlabel!("y, m")
ylabel!("Velocity in x direction, m/s")
plot!(FEM_5000.arc_length,FEM_5000.uhx, linestyle=:dot, label="FEM")
plot!(legend=:bottomright)
savefig("LidDriven_Re_5000.pdf")


#Re 10000
Bmk_10000 = DataFrame(XLSX.readtable("Benchmark/LidDriven_ux_B.xlsx", "Re10000")...)
FEM_10000 = DataFrame(XLSX.readtable("Results/Re10000_ux.xlsx", "Sheet1")...)
plot(title="Lid Driven Cavity Flow, Re=10000", Bmk_10000.y, Bmk_10000.x0p5, seriestype = :scatter, label="Benchmark")
xlabel!("y, m")
ylabel!("Velocity in x direction, m/s")
plot!(FEM_10000.arc_length,FEM_10000.uhx, linestyle=:dot, label="FEM")
plot!(legend=:bottomright)
savefig("LidDriven_Re_10000.pdf")


#SUPG with h const

#Re 1000
Bmk_1000 = DataFrame(XLSX.readtable("Benchmark/LidDriven_ux_B.xlsx", "Re1000")...)
FEM_1000h = DataFrame(XLSX.readtable("Results/Re1000_ux_hcost.xlsx", "Sheet1")...)
plot(title="Lid Driven Cavity Flow, Re=1000, SUPG h cost", Bmk_1000.y, Bmk_1000.x0p5, seriestype = :scatter, label="Benchmark")
xlabel!("y, m")
ylabel!("Velocity in x direction, m/s")
plot!(FEM_1000h.arc_length,FEM_1000h.uhx, linestyle=:dot, label="FEM")
plot!(legend=:bottomright)
savefig("LidDriven_Re_1000_hcost.pdf")

err_p = abs.(FEM_1000h.uhx - FEM_1000.uhx) ./ abs.(FEM_1000.uhx)
mean(err_p[2:end])


#Re 5000
Bmk_5000 = DataFrame(XLSX.readtable("Benchmark/LidDriven_ux_B.xlsx", "Re5000")...)
FEM_5000h = DataFrame(XLSX.readtable("Results/Re5000_ux_hcost.xlsx", "Sheet1")...)
plot(title="Lid Driven Cavity Flow, Re=5000, SUPG h cost", Bmk_5000.y, Bmk_5000.x0p5, seriestype = :scatter, label="Benchmark")
xlabel!("y, m")
ylabel!("Velocity in x direction, m/s")
plot!(FEM_5000h.arc_length,FEM_5000h.uhx, linestyle=:dot, label="FEM")
plot!(legend=:bottomright)

savefig("LidDriven_Re_5000_hcost.pdf")

err_p = abs.(FEM_5000h.uhx - FEM_5000.uhx) ./ abs.(FEM_5000.uhx)

mean(err_p[2:end])


#Re 10000
Bmk_10000 = DataFrame(XLSX.readtable("Benchmark/LidDriven_ux_B.xlsx", "Re10000")...)
FEM_10000h = DataFrame(XLSX.readtable("Results/Re10000_ux_hcost.xlsx", "Sheet1")...)
plot(title="Lid Driven Cavity Flow, Re=10000, SUPG h cost", Bmk_10000.y, Bmk_10000.x0p5, seriestype = :scatter, label="Benchmark")
xlabel!("y, m")
ylabel!("Velocity in x direction, m/s")
plot!(FEM_10000h.arc_length,FEM_10000h.uhx, linestyle=:dot, label="FEM")
plot!(legend=:bottomright)
err_p = abs.(FEM_10000h.uhx - FEM_10000.uhx) ./ abs.(FEM_10000.uhx)
mean(err_p[2:end])