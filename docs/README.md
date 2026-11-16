# QURS Shiny application — GitHub Pages

The source code lives at [Shiny_App/app.R](../Shiny_App/app.R).

The GitHub Actions workflow `.github/workflows/build-shinylive.yml` exports the Shiny application to static Shinylive assets and deploys those assets directly to GitHub Pages. The generated files are deployment artifacts and do not need to be committed to `main/docs`.

## Enable the website

1. Open **Settings → Pages** in this repository.
2. Under **Build and deployment → Source**, select **GitHub Actions** (not `Deploy from a branch`).
3. Open **Actions → Publish QURS Shinylive on GitHub Pages** and use **Run workflow** if needed.
4. Check the workflow succeeds before using the website.

Expected website URL: https://gpem-ufsm.github.io/QURS-Distribution/

The original `app.R` remains unchanged. The Monte Carlo `.rds` outputs are not needed for this interactive visualization.
