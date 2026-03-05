package main

import (
	"bytes"
	"strings"
	"testing"
)

// executeCommand runs the root command with the given args and returns stdout, stderr, and error.
func executeCommand(args ...string) (string, string, error) {
	cmd := newRootCmd()
	stdout := new(bytes.Buffer)
	stderr := new(bytes.Buffer)
	cmd.SetOut(stdout)
	cmd.SetErr(stderr)
	cmd.SetArgs(args)
	err := cmd.Execute()
	return stdout.String(), stderr.String(), err
}

func TestHelpFlag(t *testing.T) {
	out, _, err := executeCommand("--help")
	if err != nil {
		t.Fatalf("--help failed: %v", err)
	}
	if !strings.Contains(out, "vibe-vep") {
		t.Error("--help output should contain 'vibe-vep'")
	}
	if !strings.Contains(out, "annotate") {
		t.Error("--help output should list annotate subcommand")
	}
}

func TestVersionFlag(t *testing.T) {
	out, _, err := executeCommand("--version")
	if err != nil {
		t.Fatalf("--version failed: %v", err)
	}
	if !strings.Contains(out, "dev") {
		t.Errorf("--version should contain 'dev', got: %s", out)
	}
}

func TestSubcommandHelp(t *testing.T) {
	subcommands := []string{"annotate", "compare", "convert", "download", "export", "prepare", "version"}

	for _, sub := range subcommands {
		t.Run(sub, func(t *testing.T) {
			out, _, err := executeCommand(sub, "--help")
			if err != nil {
				t.Fatalf("%s --help failed: %v", sub, err)
			}
			if out == "" {
				t.Errorf("%s --help produced no output", sub)
			}
		})
	}
}

func TestAnnotateSubcommandHelp(t *testing.T) {
	for _, sub := range []string{"maf", "vcf", "variant"} {
		t.Run(sub, func(t *testing.T) {
			out, _, err := executeCommand("annotate", sub, "--help")
			if err != nil {
				t.Fatalf("annotate %s --help failed: %v", sub, err)
			}
			if !strings.Contains(out, "--assembly") {
				t.Errorf("annotate %s --help should mention --assembly flag", sub)
			}
		})
	}
}

func TestBadAssembly(t *testing.T) {
	// annotate maf with bad assembly should fail with a helpful message.
	// Use an existing file so the parser doesn't fail first.
	_, _, err := executeCommand("annotate", "maf", "--assembly", "hg99", "../../testdata/sample.maf")
	if err == nil {
		t.Fatal("expected error for bad assembly")
	}
	if !strings.Contains(err.Error(), "unsupported assembly") {
		t.Errorf("expected 'unsupported assembly' in error, got: %v", err)
	}
}

func TestMissingInputFile(t *testing.T) {
	_, _, err := executeCommand("annotate", "maf", "nonexistent_file.maf")
	if err == nil {
		t.Fatal("expected error for missing input file")
	}
}

func TestConvertMissingInput(t *testing.T) {
	// vcf2maf with missing file
	_, _, err := executeCommand("convert", "vcf2maf", "nonexistent.vcf")
	if err == nil {
		t.Fatal("expected error for missing input file")
	}
}

func TestExportParquetMissingInput(t *testing.T) {
	// export parquet without --from-cache and no input file
	_, _, err := executeCommand("export", "parquet")
	if err == nil {
		t.Fatal("expected error for missing input")
	}
	if !strings.Contains(err.Error(), "input file required") {
		t.Errorf("expected 'input file required' in error, got: %v", err)
	}
}

func TestNormalizeAssembly(t *testing.T) {
	tests := []struct {
		input   string
		want    string
		wantErr bool
	}{
		{"GRCh38", "GRCh38", false},
		{"GRCh37", "GRCh37", false},
		{"grch38", "GRCh38", false},
		{"grch37", "GRCh37", false},
		{"hg38", "GRCh38", false},
		{"hg19", "GRCh37", false},
		{"HG38", "GRCh38", false},
		{"HG19", "GRCh37", false},
		{"hg99", "", true},
		{"mm10", "", true},
		{"", "", true},
	}

	for _, tt := range tests {
		t.Run(tt.input, func(t *testing.T) {
			got, err := normalizeAssembly(tt.input)
			if (err != nil) != tt.wantErr {
				t.Errorf("normalizeAssembly(%q) error = %v, wantErr %v", tt.input, err, tt.wantErr)
				return
			}
			if got != tt.want {
				t.Errorf("normalizeAssembly(%q) = %q, want %q", tt.input, got, tt.want)
			}
		})
	}
}

func TestPickAndMostSevereMutualExclusion(t *testing.T) {
	// --pick and --most-severe together should fail
	_, _, err := executeCommand("annotate", "maf", "--pick", "--most-severe", "input.maf")
	if err == nil {
		t.Fatal("expected error for --pick + --most-severe")
	}
	if !strings.Contains(err.Error(), "mutually exclusive") {
		t.Errorf("expected 'mutually exclusive' in error, got: %v", err)
	}
}
