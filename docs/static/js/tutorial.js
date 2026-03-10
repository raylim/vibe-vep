// tutorial.js — guided tutorial sidebar for the vibe-vep interactive demo.
(function () {
  "use strict";

  const steps = [
    {
      title: "Welcome",
      content:
        "<p><strong>vibe-vep</strong> is a fast variant effect predictor " +
        "for cancer genomics. This interactive tutorial lets you try it " +
        "right in your browser.</p>" +
        "<p>The demo includes KRAS and EGFR transcript data (GRCh38).</p>" +
        "<p>Click <strong>Next</strong> to get started, or type commands " +
        "directly in the terminal.</p>",
      command: null,
    },
    {
      title: "Annotate KRAS G12C (genomic)",
      content:
        "<p>Let's annotate the most common KRAS mutation in lung cancer " +
        "using genomic coordinates.</p>" +
        "<p><code>12:25245351:C:A</code> means chromosome 12, position " +
        "25,245,351, reference C, alternate A.</p>",
      command: "vibe-vep annotate variant 12:25245351:C:A",
    },
    {
      title: "Understanding the output",
      content:
        "<p>The output shows annotations for each overlapping transcript:</p>" +
        "<ul>" +
        "<li><strong>Gene</strong> — gene symbol (KRAS)</li>" +
        "<li><strong>Transcript</strong> — Ensembl transcript ID</li>" +
        "<li><strong>Canon</strong> — YES if canonical transcript</li>" +
        "<li><strong>Consequence</strong> — SO term (e.g. missense_variant)</li>" +
        "<li><strong>Impact</strong> — HIGH, MODERATE, LOW, or MODIFIER</li>" +
        "<li><strong>HGVSc</strong> — coding DNA notation</li>" +
        "<li><strong>HGVSp</strong> — protein notation (e.g. p.Gly12Cys)</li>" +
        "</ul>" +
        "<p>Click <strong>Next</strong> to try protein notation.</p>",
      command: null,
    },
    {
      title: "Annotate with protein notation",
      content:
        "<p>You can also use protein-level notation. vibe-vep will " +
        "reverse-map it to genomic coordinates using the canonical " +
        "transcript.</p>" +
        "<p><code>KRAS G12C</code> means glycine (G) to cysteine (C) " +
        "at position 12 of the KRAS protein.</p>",
      command: "vibe-vep annotate variant KRAS G12C",
    },
    {
      title: "View a VCF file",
      content:
        "<p>The demo includes an example VCF file with variants on " +
        "chromosomes 7 and 12. Let's look at it.</p>",
      command: "cat example.vcf",
    },
    {
      title: "Annotate a VCF file",
      content:
        "<p>Now let's annotate all variants in the VCF file at once. " +
        "This processes each variant and shows the predicted consequences.</p>",
      command: "vibe-vep annotate vcf example.vcf",
    },
    {
      title: "Annotate a MAF file",
      content:
        "<p>vibe-vep also supports MAF (Mutation Annotation Format), " +
        "the standard format used by cBioPortal and TCGA. Let's view " +
        "and annotate the example MAF file.</p>",
      command: "vibe-vep annotate maf example.maf",
    },
    {
      title: "Explore on your own",
      content:
        "<p>You've completed the tutorial! Here are some things to try:</p>" +
        "<ul>" +
        '<li><code>vibe-vep genes</code> — list available genes</li>' +
        '<li><code>vibe-vep annotate variant KRAS G12D</code> — another KRAS hotspot</li>' +
        '<li><code>vibe-vep annotate variant EGFR L858R</code> — common EGFR mutation</li>' +
        '<li><code>vibe-vep annotate variant KRAS c.35G>T</code> — HGVSc notation</li>' +
        '<li><code>vibe-vep version</code> — version info</li>' +
        "</ul>" +
        "<p>See the <a href='/vibe-vep/docs/how-to-use/getting-started/'>Getting Started</a> " +
        "guide to install the full CLI.</p>",
      command: "vibe-vep genes",
    },
  ];

  let currentStep = 0;

  function init() {
    const sidebar = document.getElementById("tutorial-sidebar");
    if (!sidebar) return;

    renderStep();
  }

  function renderStep() {
    const sidebar = document.getElementById("tutorial-sidebar");
    if (!sidebar) return;

    const step = steps[currentStep];

    let html = '<div class="tutorial-step">';
    html += '<div class="tutorial-step-header">';
    html +=
      '<span class="tutorial-step-number">Step ' +
      (currentStep + 1) +
      " of " +
      steps.length +
      "</span>";
    html += "<h3>" + step.title + "</h3>";
    html += "</div>";
    html += '<div class="tutorial-step-content">' + step.content + "</div>";

    if (step.command) {
      html += '<div class="tutorial-command">';
      html +=
        '<button class="tutorial-run-btn" onclick="window.tutorialRunCommand()">';
      html += '<span class="tutorial-run-icon">&#9654;</span> Run: ';
      html += "<code>" + escapeHtml(step.command) + "</code>";
      html += "</button>";
      html += "</div>";
    }

    html += '<div class="tutorial-nav">';
    if (currentStep > 0) {
      html +=
        '<button class="tutorial-nav-btn tutorial-prev" onclick="window.tutorialPrev()">&#8592; Previous</button>';
    }
    if (currentStep < steps.length - 1) {
      html +=
        '<button class="tutorial-nav-btn tutorial-next" onclick="window.tutorialNext()">Next &#8594;</button>';
    }
    html += "</div>";
    html += "</div>";

    // Progress bar.
    html += '<div class="tutorial-progress">';
    for (let i = 0; i < steps.length; i++) {
      const cls = i === currentStep ? "active" : i < currentStep ? "done" : "";
      html += '<div class="tutorial-progress-dot ' + cls + '"></div>';
    }
    html += "</div>";

    sidebar.innerHTML = html;
  }

  function escapeHtml(text) {
    const div = document.createElement("div");
    div.textContent = text;
    return div.innerHTML;
  }

  window.tutorialRunCommand = function () {
    const step = steps[currentStep];
    if (step.command && window.vibeVEPTerminal) {
      window.vibeVEPTerminal.runCommand(step.command);
    }
  };

  window.tutorialNext = function () {
    if (currentStep < steps.length - 1) {
      currentStep++;
      renderStep();
    }
  };

  window.tutorialPrev = function () {
    if (currentStep > 0) {
      currentStep--;
      renderStep();
    }
  };

  if (document.readyState === "loading") {
    document.addEventListener("DOMContentLoaded", init);
  } else {
    init();
  }
})();
