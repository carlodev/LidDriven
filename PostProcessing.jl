using JLD2
using DataFrames, XLSX, CSV
using Plots
using Statistics

#Reading Benchmark
Bmk_1000 = DataFrame(XLSX.readtable("Benchmark/LidDriven_ux_B.xlsx", "Re1000")...)
Bmk_5000 = DataFrame(XLSX.readtable("Benchmark/LidDriven_ux_B.xlsx", "Re5000")...)
Bmk_10000 = DataFrame(XLSX.readtable("Benchmark/LidDriven_ux_B.xlsx", "Re10000")...)


#Reading SUPG Results
FEM_1000 = DataFrame(XLSX.readtable("Results/Lid_Driven_Results.xlsx", "Re_1000_SUPG")...)
FEM_5000 = DataFrame(XLSX.readtable("Results/Lid_Driven_Results.xlsx", "Re_5000_SUPG")...)
FEM_10000 = DataFrame(XLSX.readtable("Results/Lid_Driven_Results.xlsx", "Re_10000_SUPG")...)

#Reading SUPG h const Results
FEM_1000h = DataFrame(XLSX.readtable("Results/Lid_Driven_Results.xlsx", "Re_1000_SUPG_hcost")...)
FEM_5000h = DataFrame(XLSX.readtable("Results/Lid_Driven_Results.xlsx", "Re_5000_SUPG_hcost")...)
FEM_10000h = DataFrame(XLSX.readtable("Results/Lid_Driven_Results.xlsx", "Re_10000_SUPG_hcost")...)

#Reading VMS Results
VMS_1000 = DataFrame(XLSX.readtable("Results/Lid_Driven_Results.xlsx", "Re_1000_VMS")...)
VMS_5000 = DataFrame(XLSX.readtable("Results/Lid_Driven_Results.xlsx", "Re_5000_VMS")...)
VMS_10000 = DataFrame(XLSX.readtable("Results/Lid_Driven_Results.xlsx", "Re_10000_VMS")...)



#Computing Differences h const
err_p = abs.(FEM_1000h.uhx - FEM_1000.uhx) ./ abs.(FEM_1000.uhx)
mean(err_p[2:end-1])

err_p = abs.(FEM_5000h.uhx - FEM_5000.uhx) ./ abs.(FEM_5000.uhx)
mean(err_p[2:end])

err_p = abs.(FEM_10000h.uhx - FEM_10000.uhx) ./ abs.(FEM_10000.uhx)
mean(err_p[2:end-1])


#Computing Differences VMS
err_p = abs.(FEM_1000.uhx - VMS_1000.uhx) ./ abs.(FEM_1000.uhx)
mean(err_p[2:end-1])

err_p = abs.(FEM_5000.uhx - VMS_5000.uhx) ./ abs.(FEM_5000.uhx)
mean(err_p[2:end])

err_p = abs.(FEM_10000.uhx - VMS_10000.uhx) ./ abs.(FEM_10000.uhx)
mean(err_p[2:end-1])

#Plotting SUPG Results

#Re 1000
plot(title="Lid Driven Cavity Flow, Re=1000", Bmk_1000.y, Bmk_1000.x0p5, seriestype = :scatter, label="Benchmark")
xlabel!("y, m")
ylabel!("Velocity in x direction, m/s")
plot!(FEM_1000.arc_length,FEM_1000.uhx, linestyle=:dot, label="FEM")
plot!(legend=:bottomright)
savefig("LidDriven_Re_1000.pdf")


#Re 5000
plot(title="Lid Driven Cavity Flow, Re=5000", Bmk_5000.y, Bmk_5000.x0p5, seriestype = :scatter, label="Benchmark")
xlabel!("y, m")
ylabel!("Velocity in x direction, m/s")
plot!(FEM_5000.arc_length,FEM_5000.uhx, linestyle=:dot, label="FEM")
plot!(legend=:bottomright)
savefig("LidDriven_Re_5000.pdf")


#Re 10000
plot(title="Lid Driven Cavity Flow, Re=10000", Bmk_10000.y, Bmk_10000.x0p5, seriestype = :scatter, label="Benchmark")
xlabel!("y, m")
ylabel!("Velocity in x direction, m/s")
plot!(FEM_10000.arc_length,FEM_10000.uhx, linestyle=:dot, label="FEM")
plot!(legend=:bottomright)
savefig("LidDriven_Re_10000.pdf")


#Plotting SUPG with h const

#Re 1000
plot(title="Lid Driven Cavity Flow, Re=1000, SUPG h cost", Bmk_1000.y, Bmk_1000.x0p5, seriestype = :scatter, label="Benchmark")
xlabel!("y, m")
ylabel!("Velocity in x direction, m/s")
plot!(FEM_1000h.arc_length,FEM_1000h.uhx, linestyle=:dot, label="FEM")
plot!(legend=:bottomright)
savefig("LidDriven_Re_1000_hcost.pdf")

#Re 5000
plot(title="Lid Driven Cavity Flow, Re=5000, SUPG h cost", Bmk_5000.y, Bmk_5000.x0p5, seriestype = :scatter, label="Benchmark")
xlabel!("y, m")
ylabel!("Velocity in x direction, m/s")
plot!(FEM_5000h.arc_length,FEM_5000h.uhx, linestyle=:dot, label="FEM")
plot!(legend=:bottomright)
savefig("LidDriven_Re_5000_hcost.pdf")

#Re 10000
plot(title="Lid Driven Cavity Flow, Re=10000, SUPG h cost", Bmk_10000.y, Bmk_10000.x0p5, seriestype = :scatter, label="Benchmark")
xlabel!("y, m")
ylabel!("Velocity in x direction, m/s")
plot!(FEM_10000h.arc_length,FEM_10000h.uhx, linestyle=:dot, label="FEM")
plot!(legend=:bottomright)
savefig("LidDriven_Re_10000_hcost.pdf")



#Plotting VMS Results

#Re 1000
plot(title="Lid Driven Cavity Flow, Re=1000", Bmk_1000.y, Bmk_1000.x0p5, seriestype = :scatter, label="Benchmark")
xlabel!("y, m")
ylabel!("Velocity in x direction, m/s")
plot!(VMS_1000.arc_length,VMS_1000.uhx, linestyle=:dot, label="FEM")
plot!(legend=:bottomright)
savefig("LidDriven_Re_1000_VMS.pdf")


#Re 5000
plot(title="Lid Driven Cavity Flow, Re=5000", Bmk_5000.y, Bmk_5000.x0p5, seriestype = :scatter, label="Benchmark")
xlabel!("y, m")
ylabel!("Velocity in x direction, m/s")
plot!(VMS_5000.arc_length,VMS_5000.uhx, linestyle=:dot, label="FEM")
plot!(legend=:bottomright)
savefig("LidDriven_Re_5000_VMS.pdf")


#Re 10000
plot(title="Lid Driven Cavity Flow, Re=10000", Bmk_10000.y, Bmk_10000.x0p5, seriestype = :scatter, label="Benchmark")
xlabel!("y, m")
ylabel!("Velocity in x direction, m/s")
plot!(VMS_10000.arc_length,VMS_10000.uhx, linestyle=:dot, label="FEM")
plot!(legend=:bottomright)
savefig("LidDriven_Re_10000_VMS.pdf")


#Plotting VMS vs SUPG Results

#Re 1000
plot(title="Lid Driven Cavity Flow, Re=1000", Bmk_1000.y, Bmk_1000.x0p5, seriestype = :scatter, label="Benchmark")
xlabel!("y, m")
ylabel!("Velocity in x direction, m/s")
plot!(VMS_1000.arc_length,VMS_1000.uhx, linestyle=:dot, label="VMS")
plot!(FEM_1000.arc_length,FEM_1000.uhx, linestyle=:dot, label="SUPG")
plot!(legend=:bottomright)
savefig("LidDriven_Re_1000_VMSvsSUPG.pdf")


#Re 5000
plot(title="Lid Driven Cavity Flow, Re=5000", Bmk_5000.y, Bmk_5000.x0p5, seriestype = :scatter, label="Benchmark")
xlabel!("y, m")
ylabel!("Velocity in x direction, m/s")
plot!(VMS_5000.arc_length,VMS_5000.uhx, linestyle=:dot, label="VMS")
plot!(FEM_5000.arc_length,FEM_5000.uhx, linestyle=:dot, label="SUPG")
plot!(legend=:bottomright)
savefig("LidDriven_Re_5000_VMSvsSUPG.pdf")


#Re 10000
plot(title="Lid Driven Cavity Flow, Re=10000", Bmk_10000.y, Bmk_10000.x0p5, seriestype = :scatter, label="Benchmark")
xlabel!("y, m")
ylabel!("Velocity in x direction, m/s")
plot!(VMS_10000.arc_length,VMS_10000.uhx, linestyle=:dot, label="VMS")
plot!(FEM_10000.arc_length,FEM_10000.uhx, linestyle=:dot, label="SUPG")
plot!(legend=:bottomright)
savefig("LidDriven_Re_10000_VMSvsSUPG.pdf")