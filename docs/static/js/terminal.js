// terminal.js — xterm.js shell with WASM glue for vibe-vep interactive tutorial.
(function () {
  "use strict";

  // Virtual filesystem for example files.
  const virtualFS = {
    "example.vcf": null, // loaded from WASM
    "example.maf": null, // loaded from WASM
  };

  // Command history.
  const history = [];
  let historyIndex = -1;

  // Terminal state.
  let inputBuffer = "";
  let cursorPos = 0;
  let term = null;
  let wasmReady = false;
  let wasmLoading = false;

  const PROMPT = "\x1b[32mvibe-vep-demo\x1b[0m $ ";

  function initTerminal() {
    const container = document.getElementById("terminal-container");
    if (!container) return;

    term = new Terminal({
      cursorBlink: true,
      fontSize: 14,
      fontFamily: '"Fira Code", "Cascadia Code", "JetBrains Mono", Menlo, monospace',
      theme: {
        background: "#1a1b26",
        foreground: "#c0caf5",
        cursor: "#c0caf5",
        cursorAccent: "#1a1b26",
        selectionBackground: "#33467c",
        black: "#15161e",
        red: "#f7768e",
        green: "#9ece6a",
        yellow: "#e0af68",
        blue: "#7aa2f7",
        magenta: "#bb9af7",
        cyan: "#7dcfff",
        white: "#a9b1d6",
        brightBlack: "#414868",
        brightRed: "#f7768e",
        brightGreen: "#9ece6a",
        brightYellow: "#e0af68",
        brightBlue: "#7aa2f7",
        brightMagenta: "#bb9af7",
        brightCyan: "#7dcfff",
        brightWhite: "#c0caf5",
      },
      allowProposedApi: true,
    });

    const fitAddon = new FitAddon.FitAddon();
    term.loadAddon(fitAddon);
    term.open(container);
    fitAddon.fit();

    window.addEventListener("resize", () => fitAddon.fit());
    new ResizeObserver(() => fitAddon.fit()).observe(container);

    term.onData(onData);

    writeln("\x1b[1;36mWelcome to vibe-vep interactive tutorial!\x1b[0m");
    writeln("");
    loadWasm();
  }

  function loadWasm() {
    if (wasmLoading) return;
    wasmLoading = true;

    writeln("Loading annotation engine...");

    window.onVibeVEPReady = function () {
      wasmReady = true;
      // Load example VCF into virtual FS.
      if (window.vibeVEP) {
        if (window.vibeVEP.exampleVCF) {
          virtualFS["example.vcf"] = window.vibeVEP.exampleVCF();
        }
        if (window.vibeVEP.exampleMAF) {
          virtualFS["example.maf"] = window.vibeVEP.exampleMAF();
        }
      }
      writeln("\x1b[32mReady!\x1b[0m Type \x1b[1mhelp\x1b[0m to get started.\n");
      writePrompt();
    };

    const go = new Go();
    const wasmPath =
      document.getElementById("terminal-container")?.dataset.wasmPath ||
      "/vibe-vep/wasm/vibe-vep.wasm";

    WebAssembly.instantiateStreaming(fetch(wasmPath), go.importObject)
      .then((result) => {
        go.run(result.instance);
      })
      .catch((err) => {
        writeln("\x1b[31mFailed to load WASM: " + err.message + "\x1b[0m");
        writeln(
          "Make sure vibe-vep.wasm is built. See the development docs for instructions."
        );
        writePrompt();
      });
  }

  function onData(data) {
    for (let i = 0; i < data.length; i++) {
      const ch = data[i];
      const code = ch.charCodeAt(0);

      if (code === 13) {
        // Enter
        term.write("\r\n");
        const cmd = inputBuffer.trim();
        inputBuffer = "";
        cursorPos = 0;
        if (cmd) {
          history.push(cmd);
          historyIndex = history.length;
          executeCommand(cmd);
        } else {
          writePrompt();
        }
      } else if (code === 127 || code === 8) {
        // Backspace
        if (cursorPos > 0) {
          inputBuffer =
            inputBuffer.slice(0, cursorPos - 1) + inputBuffer.slice(cursorPos);
          cursorPos--;
          redrawInput();
        }
      } else if (code === 27) {
        // Escape sequences
        if (data[i + 1] === "[") {
          const seq = data[i + 2];
          if (seq === "A") {
            // Up arrow — history
            if (historyIndex > 0) {
              historyIndex--;
              inputBuffer = history[historyIndex] || "";
              cursorPos = inputBuffer.length;
              redrawInput();
            }
            i += 2;
          } else if (seq === "B") {
            // Down arrow — history
            if (historyIndex < history.length) {
              historyIndex++;
              inputBuffer = history[historyIndex] || "";
              cursorPos = inputBuffer.length;
              redrawInput();
            }
            i += 2;
          } else if (seq === "C") {
            // Right arrow
            if (cursorPos < inputBuffer.length) {
              cursorPos++;
              term.write("\x1b[C");
            }
            i += 2;
          } else if (seq === "D") {
            // Left arrow
            if (cursorPos > 0) {
              cursorPos--;
              term.write("\x1b[D");
            }
            i += 2;
          } else if (seq === "3" && data[i + 3] === "~") {
            // Delete key
            if (cursorPos < inputBuffer.length) {
              inputBuffer =
                inputBuffer.slice(0, cursorPos) +
                inputBuffer.slice(cursorPos + 1);
              redrawInput();
            }
            i += 3;
          } else {
            i += 2;
          }
        }
      } else if (code === 3) {
        // Ctrl+C
        inputBuffer = "";
        cursorPos = 0;
        term.write("^C\r\n");
        writePrompt();
      } else if (code === 12) {
        // Ctrl+L — clear
        term.clear();
        writePrompt();
        redrawInput();
      } else if (code >= 32) {
        // Printable character
        inputBuffer =
          inputBuffer.slice(0, cursorPos) + ch + inputBuffer.slice(cursorPos);
        cursorPos++;
        redrawInput();
      }
    }
  }

  function redrawInput() {
    // Clear line after prompt and rewrite.
    term.write("\r" + PROMPT + "\x1b[K" + inputBuffer);
    // Move cursor to correct position.
    const back = inputBuffer.length - cursorPos;
    if (back > 0) {
      term.write("\x1b[" + back + "D");
    }
  }

  function writePrompt() {
    term.write(PROMPT);
  }

  function writeln(text) {
    // xterm.js needs \r\n for proper line breaks; WASM output uses bare \n.
    term.write(text.replace(/\n/g, "\r\n") + "\r\n");
  }

  function executeCommand(cmd) {
    const parts = cmd.split(/\s+/);
    const base = parts[0].toLowerCase();

    if (base === "clear") {
      term.clear();
      writePrompt();
      return;
    }

    if (base === "help") {
      if (wasmReady) {
        writeln(window.vibeVEP.help());
      } else {
        writeln("WASM not loaded yet. Please wait...");
      }
      writePrompt();
      return;
    }

    if (base === "ls") {
      const files = Object.keys(virtualFS);
      writeln(files.join("  "));
      writePrompt();
      return;
    }

    if (base === "cat") {
      const file = parts[1];
      if (!file) {
        writeln("Usage: cat <file>");
      } else if (virtualFS[file] !== undefined && virtualFS[file] !== null) {
        writeln(virtualFS[file]);
      } else {
        writeln("cat: " + file + ": No such file");
      }
      writePrompt();
      return;
    }

    if (base === "vibe-vep") {
      if (!wasmReady) {
        writeln("WASM not loaded yet. Please wait...");
        writePrompt();
        return;
      }

      const sub = (parts[1] || "").toLowerCase();

      if (sub === "version") {
        writeln(window.vibeVEP.version());
        writePrompt();
        return;
      }

      if (sub === "help" || !sub) {
        writeln(window.vibeVEP.help());
        writePrompt();
        return;
      }

      if (sub === "genes") {
        writeln(window.vibeVEP.listGenes());
        writePrompt();
        return;
      }

      if (sub === "annotate") {
        const subSub = (parts[2] || "").toLowerCase();

        if (subSub === "variant") {
          const spec = parts.slice(3).join(" ");
          if (!spec) {
            writeln("Usage: vibe-vep annotate variant <spec>");
            writeln("  e.g. vibe-vep annotate variant 12:25245351:C:A");
            writeln("  e.g. vibe-vep annotate variant KRAS G12C");
          } else {
            writeln(window.vibeVEP.annotateVariant(spec));
          }
          writePrompt();
          return;
        }

        if (subSub === "vcf") {
          const file = parts[3];
          if (!file) {
            writeln("Usage: vibe-vep annotate vcf <file>");
          } else if (
            virtualFS[file] !== undefined &&
            virtualFS[file] !== null
          ) {
            writeln(window.vibeVEP.annotateVCF(virtualFS[file]));
          } else {
            writeln("Error: file not found: " + file);
          }
          writePrompt();
          return;
        }

        if (subSub === "maf") {
          const file = parts[3];
          if (!file) {
            writeln("Usage: vibe-vep annotate maf <file>");
          } else if (
            virtualFS[file] !== undefined &&
            virtualFS[file] !== null
          ) {
            writeln(window.vibeVEP.annotateMAF(virtualFS[file]));
          } else {
            writeln("Error: file not found: " + file);
          }
          writePrompt();
          return;
        }

        writeln("Usage: vibe-vep annotate variant|vcf|maf ...");
        writePrompt();
        return;
      }

      writeln("Unknown subcommand: " + sub);
      writeln('Type "vibe-vep help" for usage.');
      writePrompt();
      return;
    }

    writeln("Command not found: " + base);
    writeln('Type "help" for available commands.');
    writePrompt();
  }

  // Public API for tutorial.js to inject commands.
  window.vibeVEPTerminal = {
    runCommand: function (cmd) {
      if (!term) return;
      inputBuffer = cmd;
      cursorPos = cmd.length;
      redrawInput();
      // Simulate Enter.
      term.write("\r\n");
      inputBuffer = "";
      cursorPos = 0;
      history.push(cmd);
      historyIndex = history.length;
      executeCommand(cmd);
    },
    isReady: function () {
      return wasmReady;
    },
  };

  // Initialize when DOM is ready.
  if (document.readyState === "loading") {
    document.addEventListener("DOMContentLoaded", initTerminal);
  } else {
    initTerminal();
  }
})();
