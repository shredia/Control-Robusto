%% ========================================================================
% LOGGEAR SOLO SEÑALES QUE TERMINAN EN _M1 O _M2
% ========================================================================

modelo = 'Simulacion_stepper';

load_system(modelo);


%% ========================================================================
% LIMPIAR INSTRUMENTACIÓN ANTERIOR
% ========================================================================

set_param(modelo,'InstrumentedSignals',[]);


%% ========================================================================
% ACTIVAR SIGNAL LOGGING GLOBAL
% ========================================================================

set_param(modelo, ...
    'SignalLogging','on', ...
    'SignalLoggingName','logsout');


%% ========================================================================
% BUSCAR TODAS LAS LÍNEAS
% ========================================================================

lineas = find_system( ...
    modelo, ...
    'FollowLinks','on', ...
    'LookUnderMasks','all', ...
    'FindAll','on', ...
    'Type','line');


contador = 0;
nombres_loggeados = strings(0);


fprintf('\n');
fprintf('=============================================================\n');
fprintf(' ACTIVANDO SIGNAL LOGGING\n');
fprintf(' SOLO SEÑALES *_M1 Y *_M2\n');
fprintf('=============================================================\n');


%% ========================================================================
% RECORRER TODAS LAS LÍNEAS
% ========================================================================

for k = 1:numel(lineas)

    hLine = lineas(k);


    %% --------------------------------------------------------------------
    % OBTENER NOMBRE
    % ---------------------------------------------------------------------

    try
        nombre = strtrim(get_param(hLine,'Name'));
    catch
        nombre = '';
    end


    if isempty(nombre)
        continue;
    end


    %% --------------------------------------------------------------------
    % SOLO *_M1 Y *_M2
    % ---------------------------------------------------------------------

    if ~(endsWith(nombre,'_M1') || endsWith(nombre,'_M2'))
        continue;
    end


    %% --------------------------------------------------------------------
    % EVITAR NOMBRES DUPLICADOS
    % ---------------------------------------------------------------------

    if any(strcmp(nombres_loggeados,nombre))
        continue;
    end


    %% --------------------------------------------------------------------
    % OBTENER PUERTO DE ORIGEN
    % ---------------------------------------------------------------------

    try

        hPort = get_param(hLine,'SrcPortHandle');

    catch

        hPort = -1;

    end


    if isempty(hPort) || hPort == -1
        continue;
    end


    %% --------------------------------------------------------------------
    % ACTIVAR DATA LOGGING EN EL PUERTO
    % ---------------------------------------------------------------------

    try

        set_param( ...
            hPort, ...
            'DataLogging', ...
            'on');


        contador = contador + 1;

        nombres_loggeados(end+1) = string(nombre);


        fprintf( ...
            '%3d | %s\n', ...
            contador, ...
            nombre);


    catch ME

        fprintf( ...
            'NO LOG | %-30s | %s\n', ...
            nombre, ...
            ME.message);

    end

end


%% ========================================================================
% RESUMEN
% ========================================================================

fprintf('\n');
fprintf('=============================================================\n');
fprintf(' Señales activadas: %d\n',contador);
fprintf('=============================================================\n');


save_system(modelo);