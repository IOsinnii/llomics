// Runs Python in a module Web Worker so the page stays responsive. Loads Pyodide from the site's own
// pyodide/ folder (tools/build_site.py puts a pinned release there), copies the llomics files listed
// in repo/manifest.json into Pyodide's file system with the repository's layout, imports
// web/bridge.py, and answers {id, fn, args} messages with bridge.call(fn, args).

import { loadPyodide } from "./pyodide/pyodide.mjs";

function status(text) {
  postMessage({ type: "status", text });
}

const ready = (async () => {
  status("Loading Python for the browser (about 13 MB, cached after the first visit)");
  const pyodide = await loadPyodide({ indexURL: new URL("./pyodide/", import.meta.url).href });

  status("Loading llomics");
  const manifest = await (await fetch("repo/manifest.json", { cache: "no-cache" })).json();
  await Promise.all(
    manifest.files.map(async (path) => {
      const response = await fetch("repo/" + path, { cache: "no-cache" });
      if (!response.ok) throw new Error(`could not load repo/${path} (${response.status})`);
      const target = "/repo/" + path;
      pyodide.FS.mkdirTree(target.slice(0, target.lastIndexOf("/")));
      pyodide.FS.writeFile(target, new Uint8Array(await response.arrayBuffer()));
    }),
  );
  pyodide.runPython('import sys; sys.path.insert(0, "/repo/web"); import bridge');
  const bridge = pyodide.pyimport("bridge");
  postMessage({ type: "ready", manifest });
  return bridge;
})();

ready.catch((error) => postMessage({ type: "failed", text: String(error) }));

onmessage = async (event) => {
  const { id, fn, args } = event.data;
  try {
    const bridge = await ready;
    postMessage({ id, out: bridge.call(fn, JSON.stringify(args || {})) });
  } catch (error) {
    postMessage({ id, error: String(error) });
  }
};
