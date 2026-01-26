# import Pkg; Pkg.add("Plots") # 如果还没有安装 Plots 包，请取消注释此行运行一次
using Plots
using Plots.Measures

# 1. 物理与数值参数
L = 1000.0          # 管长 (m)
D = 0.5             # 管径 (m)
a = 1000.0          # 波速 (m/s)
g = 9.81            
f = 0.02            # 摩阻系数
V0 = 2.0            # 初始流速 (m/s)
H_res = 100.0       # 水库水头 (m)
T_close = 0.01      # 快速关阀时间 (s) - 设短一点可以看到更陡峭的波

nx = 1000            
dx = L / nx
dt = dx / a         
total_time = 16.0    # 模拟时间（s）。模拟 4 秒即足以观察一个完整的往复
nt = Int(total_time ÷ dt)

observe_x = 500      # 观察点位置，取1000m时为阀门处

# 2. 初始化数组
x = 0:dx:L
V = fill(V0, nx + 1)
H = [H_res - (i-1)*dx * (f*V0^2/(2*g*D)) for i in 1:nx+1]
H0_valve = H[end]

V_new, H_new = copy(V), copy(H)
t_hist, H_valve_hist, V_valve_hist = Float64[], Float64[], Float64[]

# 3. 动画生成
anim = Animation()
for t_step in 1:nt

    # 进度显示
    if t_step % (nt ÷ 20) == 0
        @info "当前计算进度: $(round(t_step/nt*100, digits=1))%"
    end

    curr_t = t_step * dt
    
    # --- 核心计算逻辑 (MOC) ---
    for i in 2:nx
        CP = V[i-1] + (g/a)*H[i-1] - (f*dt/(2*D))*V[i-1]*abs(V[i-1])
        CM = V[i+1] - (g/a)*H[i+1] - (f*dt/(2*D))*V[i+1]*abs(V[i+1])
        H_new[i] = (CP - CM) / (2g/a)
        V_new[i] = (CP + CM) / 2
    end

    # 边界条件
    CM_up = V[2] - (g/a)*H[2] - (f*dt/(2*D))*V[2]*abs(V[2])
    H_new[1], V_new[1] = H_res, CM_up + (g/a)*H_res

    CP_down = V[nx] + (g/a)*H[nx] - (f*dt/(2*D))*V[nx]*abs(V[nx])
    tau = max(0.0, 1.0 - curr_t / T_close)
    if tau > 0
        K = (tau * V0)^2 / H0_valve
        ca = g / (a * K)
        V_new[nx+1] = (-1.0 + sqrt(1.0 + 4.0 * ca * CP_down)) / (2.0 * ca)
        H_new[nx+1] = (CP_down - V_new[nx+1]) / (g/a)
    else
        V_new[nx+1], H_new[nx+1] = 0.0, CP_down / (g/a)
    end

    V .= V_new; H .= H_new
    
    # 记录时间序列
    push!(t_hist, curr_t); push!(H_valve_hist, H[Int(observe_x/dx)+1]); push!(V_valve_hist, V[Int(observe_x/dx)+1])

    # 水头图和流速变化
    p1 = plot(x, H, title="Head Profile (t = $(round(curr_t, digits=2)))", ylabel="H (m)", color=:blue, ylims=(H_res-250, H_res+250), label="")
    p2 = plot(x, V, title="Velocity Profile (t = $(round(curr_t, digits=2)))", ylabel="V (m/s)", color=:red, ylims=(-V0-1, V0+1), label="")
    
    # observe_x处随时间的变化
    p3 = plot(t_hist, H_valve_hist, 
              ylabel="H (m)", 
              label="Head", 
              color=:blue, 
              title="History (x=$observe_x)", 
              xlabel="Time (s)", 
              ylims=(H_res-250, H_res+250), 
              legend=:topleft)
    p3 = plot!(twinx(), t_hist, V_valve_hist, 
               ylabel="V (m/s)", 
               label="Velocity", 
               color=:red, 
               ylims=(-V0-1, V0+1), 
               legend=:topright,
               grid=false)

    # 最终组合
    # 每 100 步存一帧，提升生成速度
    if t_step % 100 == 0 || t_step == nt
        p_final = plot(p1, p2, p3, layout=(3, 1), size=(1200, 800), 
                        left_margin=10mm, right_margin=20mm, top_margin=15mm, bottom_margin=10mm)
        frame(anim, p_final)
    end
end

gif(anim, "water_hammer.gif", fps = 10)