# QURS Shiny application (GitHub Pages)

This directory is the publication target for the [QURS Shiny application](../Shiny_App/app.R).

The GitHub Actions workflow `.github/workflows/build-shinylive.yml` exports the app with the R package `shinylive`, then commits the static web output here.

After the workflow has completed successfully, enable GitHub Pages in **Settings → Pages** using **Deploy from a branch**, branch `main`, folder `/docs`.

Expected URL: https://gpem-ufsm.github.io/QURS-Distribution/

The original application is retained at `Shiny_App/app.R`. No simulation output files are required to run this interactive visualization.
