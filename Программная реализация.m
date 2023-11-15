%%
clear, clc;

[l, u_max, m_0, M, step, T, eps, n] = ParamsInput();

ts = linspace(T/(n+1), n*T/(n+1), n); %Рассматриваемые моменты переключений
tspan = [0 : step : T];
if tspan(end) ~= T
    tspan(end + 1) = T;
end
 
count = size(tspan, 2);
x_memory = zeros(count - 2 , 1);
psi_memory = zeros(count - 2, 1);
u_memory = zeros(count - 2, 1);

g = 9.8;  %Гравитационная постоянная
i_best = 0;
h_best = 0;

for i = 1 : n
    tau_1 = ts(i);

    [f1, cond] = FirstStageSearchPsi_0(u_max, g, l, m_0, tau_1);

    for j = i + 1 : n
        tau_2 = ts(j);

        %Строится траектория для случая u_max -> u* -> ..., с переключениями в моменты tau_1 и tau_2
        [x_memory, psi_memory, u_memory, i_best, h_best] = trajectory('u_max -> u* -> ...', f1, cond, g, l, m_0, M, tspan, step, tau_1, tau_2, T, x_memory, psi_memory, u_memory, eps, h_best, i_best, u_max);                    

        %Строится траектория для случая u_max -> 0 -> ..., с переключениями в моменты tau_1 и tau_2
        [x_memory, psi_memory, u_memory, i_best, h_best] = trajectory('u_max -> 0 -> ...', f1, cond, g, l, m_0, M, tspan, step, tau_1, tau_2, T, x_memory, psi_memory, u_memory, eps, h_best, i_best, u_max);                    
    end

    %Случай соответствует u_max -> 0 с выходом на границу конечного мн-ва по x_2 (v(T) \in [-eps, eps]) 
    [x_memory, psi_memory, u_memory, i_best, h_best] = trajectory('u_max -> 0, eps', f1, cond, g, l, m_0, M, tspan, step, tau_1, T, T, x_memory, psi_memory, u_memory, eps, h_best, i_best, u_max);
    
    %Случай соответствует u_max -> 0 с выходом на границу конечного мн-ва по x_3 (m(T) = M) 
    [x_memory, psi_memory, u_memory, i_best, h_best] = trajectory('u_max -> 0, M', f1, cond, g, l, m_0, M, tspan, step, tau_1, T, T, x_memory, psi_memory, u_memory, eps, h_best, i_best, u_max);
end

if i_best
    x_best = x_memory( : , i_best * 3 - 1 : i_best * 3 + 1); 
    psi_best = psi_memory( : , i_best * 2 : i_best * 2 + 1);
    u_best = u_memory( : , i_best + 1);
 
    %Полученное значение функционала на оптимальной траектории
    disp(['J(u) = ', num2str(h_best)]);
else
    disp('Задача не разрешима');
end

%Визуализация результатов
draw_graphics(x_memory, psi_memory, u_memory, tspan(1 : end - 2), u_max, T, i_best);

%Сократить длину строк
%%
function [l, u_max, m_0, M, step, T, eps, n] = ParamsInput()

    l = input('Введите l > 0 - коэффициент, определяющий силу, действующую на ракету со стороны сгорающего топлива\n');
    u_max = input('Введите максимально возможную скорость подачи топлива u_max\n');
    m_0 = input('Введите начальную массу ракеты с топливом m_0\n');
    M = input('Введите массу M корпуса ракеты без топлива\n');
    step = input('Введите шаг интегрирования систем дифференциальных уравнений step\n');
    T = input('Введите конечный момент времени T\n');
    eps = input('Введите eps - ограничение на конечную скорость (v(T) \in [-eps, eps])\n');
    n = input('Введите число элементов в сетке по времени\n');
    
end

function [f1, cond] = FirstStageSearchPsi_0(u_max, g, l, m_0, tau_1)
    syms psi_2_0 psi_3_0 x_2(t) x_3(t) psi_2(t) psi_3(t) u(t) ...
         v(t) m(t) f1(psi_2_0, psi_3_0);

    u(t) = u_max;
    eqn1_uMax = diff(x_2, t) == -g + (x_2 + l) / (x_3) * u;
    eqn2_uMax = diff(x_3, t) == -u;
    eqn3_uMax = diff(psi_2, t) == -1 - psi_2 * u / x_3;
    eqn4_uMax = diff(psi_3, t) == psi_2 * (x_2 + l) / (x_3)^2 * u; 
    
    %Первый участок на u_max во всех случаях
    cond = [0, m_0, psi_2_0, psi_3_0];
    res = dsolve([eqn1_uMax, eqn2_uMax, eqn3_uMax, eqn4_uMax], ...
                 [x_2(0), x_3(0), psi_2(0), psi_3(0)] == cond);

    v(t) = res.x_2;
    m(t) = res.x_3;
    psi(t) = [res.psi_2, res.psi_3];
    f1(psi_2_0, psi_3_0) = dot(psi(tau_1), [(v(tau_1) + l) / m(tau_1), -1]);  

    cond = [v(tau_1), m(tau_1), psi(tau_1)];
end

function [x_memory, psi_memory, u_memory, i_best, h_best] = trajectory(str, f1, cond, g, l, m_0, M, tspan, step, tau_1, tau_2, T, x_memory, psi_memory, u_memory, eps, h_best, i_best, u_max)
    
    psi_0 = SecondStageSearchPsi_0(str, f1, cond, g, l, tau_1, tau_2);
    [x_memory, psi_memory, u_memory, i_best, h_best] = Check_and_Building(str, g, m_0, M, psi_0, tspan, tau_1, step, tau_2, x_memory, psi_memory, u_memory, eps, h_best, i_best, l, T, u_max);               
    
end

function psi_0 = SecondStageSearchPsi_0(str, f1, cond, g, l, tau_1, tau_2)
    %Подготовка к символьному решению дифференциальных уравнений              
    syms psi_2_0 psi_3_0 x_2(t) x_3(t) psi_2(t) psi_3(t) u(t) ...
         v(t) m(t) f2(psi_2_0, psi_3_0);

    u(t) = g * x_3(t) / (x_2(t) + l);
    eqn1_u = diff(x_2, t) == 0;
    eqn2_u = diff(x_3, t) == -u;
    eqn3_u = diff(psi_2, t) == -1 - psi_2 * u / x_3;
    eqn4_u = diff(psi_3, t) == psi_2 * (x_2 + l) / (x_3)^2 * u;

    eqn1_0 = diff(x_2, t) == -g;
    eqn2_0 = diff(x_3, t) == 0;
    eqn3_0 = diff(psi_2, t) == -1;
    eqn4_0 = diff(psi_3, t) == 0;    
    
    switch str(10)
        case 'u'
            %Второй участок - особый режим
            res = dsolve([eqn1_u, eqn2_u, eqn3_u, eqn4_u], ...
                         [x_2(tau_1), x_3(tau_1), psi_2(tau_1), psi_3(tau_1)] == cond); 

            v(t) = res.x_2;
            m(t) = res.x_3;
            psi(t) = [res.psi_2, res.psi_3];
            f2(psi_2_0, psi_3_0) = dot(psi(tau_2), [(v(tau_2) + l) / m(tau_2), -1]);  
        case '0'
            %Второй участок 0-е управление
            res = dsolve([eqn1_0, eqn2_0, eqn3_0, eqn4_0], ...
                         [x_2(tau_1), x_3(tau_1), psi_2(tau_1), psi_3(tau_1)] == cond); 

            v(t) = res.x_2;
            m(t) = res.x_3;
            psi(t) = [res.psi_2, res.psi_3];
            
            if all(str(end - 2 : end) == '...')
                %Если после 0-го идёт другой режим
                f2(psi_2_0, psi_3_0) = dot(psi(tau_2), [(v(tau_2) + l) / m(tau_2), -1]); 
            else 
                if all(str(end - 2 : end) == 'eps')
                    %Если на нулевом всё заканчивается, при этом предполагается выход на границу по x_2 в момент T
                    f2(psi_2_0, psi_3_0) = dot(psi(tau_2), [0, 1]);
                else
                    %Если на нулевом всё заканчивается, при этом предполагается выход на границу по x_3 в момент T
                    f2(psi_2_0, psi_3_0) = dot(psi(tau_2), [1, 0]);
                end
            end
    end

    sol = solve([f1(psi_2_0, psi_3_0) == 0, f2(psi_2_0, psi_3_0) == 0], [psi_2_0, psi_3_0]);
    psi_0 = double([sol.psi_2_0, sol.psi_3_0]);        
end

function [x_memory, psi_memory, u_memory, i_best, h_best] = Check_and_Building(str, g, m_0, M, psi_0, tspan, tau_1, step, tau_2, x_memory, psi_memory, u_memory, eps, h_best, i_best, l, T, u_max)
    %Системы для численного решения
    dxdt_uMax = @(t, x) [x(2);
                        -g + (x(2) + l) / x(3) * u_max;
                        -u_max;
                        -1 - x(4) * u_max / x(3);
                        x(4) * (x(2) + l) / x(3)^2 * u_max;];

    dxdt_u = @(t, x) [x(2);
                      0;
                      -g * x(3) / (x(2) + l);
                      -1 - x(4) * g / (x(2) + l);
                      x(4) / x(3) * g;]; 

    dxdt_0 = @(t, x) [x(2);
                      -g;
                      0;
                      -1;
                      0;];    
 
    delta = 0.00001;
    
    %Численно определяем траекторию для первых двух режимов исходя из найденного psi_0 и, 
    %если полученная траектория удовлетворяет заявленным режимам, то достраиваем её
    cond = [0, 0, m_0, psi_0(1), psi_0(2)];
    tspan_curr = [tspan(tspan < tau_1), tau_1];
    [~, x_first] = ode45(dxdt_uMax, tspan_curr, cond);
    
    cond = x_first(end, :);
    if tspan_curr(end) ~= step + tspan_curr(end-1)
        x_first = x_first(1 : end - 1, :);
    end
    tspan_curr = [tau_1, tspan(tspan < tau_2 & tspan > tau_1), tau_2];
    u_curr = ones(size(x_first, 1), 1) * u_max;
 
    switch str(10)
        case 'u'
            [~, x_second] = ode45(dxdt_u, tspan_curr, cond);
            u_curr = [u_curr; x_second(2 : end, 3) * g ./ (x_second(2 : end, 2) + l)];
            
            if all((x_first(1:end-1, 2) + l) ./ x_first(1:end-1, 3) .* x_first(1:end-1, 4) - x_first(1:end-1, 5) > 0) && ...
               all((x_second(2:end-1, 2) + l) ./ x_second(2:end-1, 3) .* x_second(2:end-1, 4) - x_second(2:end-1, 5) <= delta) && ...
               all((x_second(2:end-1, 2) + l) ./ x_second(2:end-1, 3) .* x_second(2:end-1, 4) - x_second(2:end-1, 5) >= -delta)
            
               %Достраиваем траекторию и сохраняем её
               [x, u] = buildTrajectory(x_first, x_second, tau_2, tspan, tspan_curr, step, T, g, l, u_max, dxdt_uMax, dxdt_u, dxdt_0, u_curr);
               x_memory = [x_memory, x(1 : size(x_memory, 1), 1 : 3)];
               psi_memory = [psi_memory, x(1 : size(psi_memory, 1), 4 : 5)]; 
               u_memory = [u_memory, u(1 : size(u_memory, 1))]; 
               
               %Запоминаем характеристики траектории, если её показатели оказались наилучшими на текущий момент
               [i_best, h_best] = best(x(end, :), h_best, i_best, eps, size(x_memory, 2), M);
            end
        case '0'
            [~, x_second] = ode45(dxdt_0, tspan_curr, cond);
            u_curr = [u_curr; x_second(2 : end, 3) * 0];
            
            if all((x_first(1:end-1, 2) + l) ./ x_first(1:end-1, 3) .* x_first(1:end-1, 4) - x_first(1:end-1, 5) > 0) && ...
               all((x_second(2:end-1, 2) + l) ./ x_second(2:end-1, 3) .* x_second(2:end-1, 4) - x_second(2:end-1, 5) < 0)

               %Достраиваем траекторию и сохраняем её
               [x, u] = buildTrajectory(x_first, x_second, tau_2, tspan, tspan_curr, step, T, g, l, u_max, dxdt_uMax, dxdt_u, dxdt_0, u_curr);
               x_memory = [x_memory, x(1 : size(x_memory, 1), 1 : 3)];
               psi_memory = [psi_memory, x(1 : size(psi_memory, 1), 4 : 5)]; 
               u_memory = [u_memory, u(1 : size(u_memory, 1))]; 

               %Запоминаем характеристики траектории, если её показатели оказались наилучшими на текущий момент
               [i_best, h_best] = best(x(end, :), h_best, i_best, eps, size(x_memory, 2), M);
            end
    end
end

function [x, u_curr] = buildTrajectory(x_first, x_second, tau_2, tspan, tspan_curr, step, T, g, l, u_max, dxdt_uMax, dxdt_u, dxdt_0, u_curr)

        delta = 0.00001; %Погрешность определения режима управления
        x = [x_first; x_second(2 : end, :)];                             
        t_curr = tau_2;

        while t_curr <= T - step
            %sgn - величина, по которой из условия максимума определяем оптимальное в данный момент управление
            sgn = (x(end, 2) + l) / x(end, 3) * x(end, 4) - x(end, 5);
            
            x0 = x(end, :);
            if tspan_curr(end) ~= step + tspan_curr(end-1) 
                x = x(1 : end - 1, :);
                u_curr = u_curr(1 : end - 1, :);
            end
            tspan_curr = [t_curr, tspan(tspan > t_curr & tspan <= T)];

            %По sgn определяем следующий режим и интегрируем траекторию в нём
            if sgn >= delta

                [tspan_curr, x_curr] = ode45(dxdt_uMax, tspan_curr, x0, odeset('Events', @(t, x) UChange(t, x, l)));
                u_curr = [u_curr; ones(size(tspan_curr, 1) - 1, 1) * u_max];
            else
                if sgn <= -delta
                
                    [tspan_curr, x_curr] = ode45(dxdt_0, tspan_curr, x0, odeset('Events', @(t, x) UChange(t, x, l)));             
                    u_curr = [u_curr; ones(size(tspan_curr, 1) - 1, 1) * 0];
                else
                    
                    [tspan_curr, x_curr] = ode45(dxdt_u, tspan_curr, x0, odeset('Events', @(t, x) SpecialOut(t, x, l, delta)));
                    u_curr = [u_curr; g * x_curr(2 : end, 3) ./ (x_curr(2 : end, 2) + l)];
                end
            end

            t_curr = tspan_curr(end);
            x = [x; x_curr(2 : end, :)];        
        end
end

function [i_best, h_best] = best(vec, h_best, i_best, eps, sz, M)
    %Проверяем попадание траектории в конечное мн-во и достижение ею наибольшей высоты
    if (vec(2) >= -eps) && (vec(2) <= eps) && ...
       (vec(1) > h_best) && (vec(3) >= M)

        i_best = (sz - 1) / 3;
        h_best = vec(1);
    end
end

%Функция события: выход из особого режима
function [value, isterminal, direction] = SpecialOut(t, x, l, eps)
    value = abs(x(4) * (x(2) + l) / x(3) - x(5)) < eps;
    isterminal = 1;   
    direction = 0;   
end          

%Функция события: смена управления
function [value, isterminal, direction] = UChange(t, x, l)
    value = x(4) * (x(2) + l) / x(3) - x(5);  
    isterminal = 1;   
    direction = 0;   
end   

function draw_graphics(x_memory, psi_memory, u_memory, t, u_max, T, i_best)

    if i_best         
        disp('Введите номера осей, графики результатов в которых вы хотите увидеть. Вторым параметром после номера для каждой оси обозначьте полноту вывода (a - только оптимальную, b -все ).')
        disp('1) (t, x_2)')
        disp('2) (t, x_3)')
        disp('3) (t, psi_2)')
        disp('4) (t, psi_3)')
        disp('5) (t, u)')
        disp('Пример ввода:')
        disp("[1 'a'; 1 'b'; 2 'a'; 2 'b'; 3 'a'; 3 'b'; 4 'a'; 4 'b'; 5 'a'; 5 'b']")
        params = input('');  
    else       
        params = [1 'b'; 2 'b'; 3 'b'; 4 'b'; 5 'b'];
    end
    
    a = params(:, 1);
    b = params(:, 2);
    F = figure;
    G = uitabgroup('Parent', F); 
    
    for i = 1 : size(params, 1)
        tab = uitab(G);
        ax = axes('parent', tab); 
        
        switch a(i)
            case 1
                name1 = 't';
                name2 = 'x_2';
                hold(tab.Children(end), 'on');
                xlim([-1, T + 1]);
                
                if (b(i) == 'a')
                    plot(tab.Children(end), t, x_memory(:, i_best * 3), 'Color', 'red', 'DisplayName', 'Opt');
                    legend(tab.Children(end), 'Opt');
                    hold(tab.Children(end), 'off');
                else
                    plot(tab.Children(end), t, x_memory(:, 3 : 3 : end), 'Color', 'blue', 'DisplayName', 'All');
                    legend(tab.Children(end), 'All');
                    if i_best
                        plot(tab.Children(end), t, x_memory(:, i_best * 3), 'Color', 'red', 'DisplayName', 'Opt');
                    end
                    hold(tab.Children(end), 'off');
                end
            case 2
                name1 = 't';
                name2 = 'x_3';
                hold(tab.Children(end), 'on');
                xlim([-1, T + 1]);
                
                if (b(i) == 'a')
                    plot(tab.Children(end), t, x_memory(:, i_best * 3 + 1), 'Color', 'red', 'DisplayName', 'Opt');
                    legend(tab.Children(end), 'Opt');
                    hold(tab.Children(end), 'off');
                else
                    plot(tab.Children(end), t, x_memory(:, 4 : 3 : end), 'Color', 'blue', 'DisplayName', 'All');
                    legend(tab.Children(end), 'All');
                    if i_best
                        plot(tab.Children(end), t, x_memory(:, i_best * 3 + 1), 'Color', 'red', 'DisplayName', 'Opt');
                    end
                    hold(tab.Children(end), 'off');
                end
            case 3
                name1 = 't';
                name2 = 'psi_2';
                hold(tab.Children(end), 'on');
                xlim([-1, T + 1]);
                
                if (b(i) == 'a')
                    plot(tab.Children(end), t, psi_memory(:, i_best * 2), 'Color', 'red', 'DisplayName', 'Opt');
                    legend(tab.Children(end), 'Opt');
                    hold(tab.Children(end), 'off');
                else
                    plot(tab.Children(end), t, psi_memory(:, 2 : 2 : end), 'Color', 'blue', 'DisplayName', 'All');
                    legend(tab.Children(end), 'All');
                    if i_best
                        plot(tab.Children(end), t, psi_memory(:, i_best * 2), 'Color', 'red', 'DisplayName', 'Opt');
                    end
                    hold(tab.Children(end), 'off');
                end
            case 4
                name1 = 't';
                name2 = 'psi_3';
                hold(tab.Children(end), 'on');
                xlim([-1, T + 1]);
                
                if (b(i) == 'a')
                    plot(tab.Children(end), t, psi_memory(:, i_best * 2 + 1), 'Color', 'red', 'DisplayName', 'Opt');
                    legend(tab.Children(end), 'Opt');
                    hold(tab.Children(end), 'off');
                else
                    plot(tab.Children(end), t, psi_memory(:, 3 : 2 : end), 'Color', 'blue', 'DisplayName', 'All');
                    legend(tab.Children(end), 'All');
                    if i_best
                        plot(tab.Children(end), t, psi_memory(:, i_best * 2 + 1), 'Color', 'red', 'DisplayName', 'Opt');
                    end
                    hold(tab.Children(end), 'off');
                end
            case 5
                name1 = 't';
                name2 = 'u';
                hold(tab.Children(end), 'on');
                xlim([-1, T + 1]);
                ylim([-1, u_max+1]);
                
                if (b(i) == 'a')
                    plot(tab.Children(end), t, u_memory(:, i_best + 1), 'Color', 'red', 'DisplayName', 'Opt');
                    legend(tab.Children(end), 'Opt');
                    hold(tab.Children(end), 'off');
                else
                    plot(tab.Children(end), t, u_memory(:, 2 : end), 'Color', 'blue', 'DisplayName', 'All');
                    legend(tab.Children(end), 'All');
                    if i_best
                        plot(tab.Children(end), t, u_memory(:, i_best + 1), 'Color', 'red', 'DisplayName', 'Opt');
                    end
                    hold(tab.Children(end), 'off');
                end
        end
        
        xlabel(tab.Children(end), name1);
        ylabel(tab.Children(end), name2); 
        
        if (b(i) == 'a')
            type = 'opt';
        else
            type = 'all';
        end      
        
        tab.Title = strcat(type, ' (', name1, ',', name2, ')');
    end
end
