%% Replace all interpreters with set value
% set_interpreters(value)
% "latex", "none", or "remove" (resets to default, which is "none");
function set_interpreters(value)
    factories = fieldnames(get(groot,"factory"));
    factories = string(factories);
    interpreters = factories(contains(factories, "Interpreter"));
    defaults = strrep(interpreters, "factory", "default");
    for terp = defaults'
            set(groot, terp, value);
    end
end
