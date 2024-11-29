omen = openrocket("data/OMEN.ork");
sim = omen.sims(1);
data = openrocket.get_data(sim, ["Altitude", "Vertical velocity", "Total velocity"]);

figure;
tiledlayout(2,1);
nexttile;
plot(data.Time, data.Altitude);
nexttile;
plot(data.Time, data.("Vertical velocity"));
plot(data.Time, data.("Total velocity"));
