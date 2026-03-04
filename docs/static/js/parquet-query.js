/**
 * VibeVEPClient — query pre-computed variant annotations from a remote Parquet
 * file using DuckDB-WASM. DuckDB makes HTTP range requests, fetching only the
 * relevant row groups (~20-50 KB per query).
 *
 * Usage:
 *   const client = new VibeVEPClient("https://example.com/annotations.parquet");
 *   await client.init();
 *   const results = await client.lookupVariant("12", 25245350, "C", "A");
 *   await client.close();
 */

const DUCKDB_CDN = "https://cdn.jsdelivr.net/npm/@duckdb/duckdb-wasm@1.29.0/dist/";

export class VibeVEPClient {
  /**
   * @param {string} parquetUrl - URL to the Parquet file (must support CORS + range requests)
   */
  constructor(parquetUrl) {
    this._parquetUrl = parquetUrl;
    this._db = null;
    this._conn = null;
  }

  /** Load DuckDB-WASM from CDN and create a view over the remote Parquet file. */
  async init() {
    // Dynamically import DuckDB-WASM
    const duckdb = await import(DUCKDB_CDN + "duckdb-esm.js");

    const bundle = await duckdb.selectBundle({
      mvp: {
        mainModule: DUCKDB_CDN + "duckdb-mvp.wasm",
        mainWorker: DUCKDB_CDN + "duckdb-browser-mvp.worker.js",
      },
      eh: {
        mainModule: DUCKDB_CDN + "duckdb-eh.wasm",
        mainWorker: DUCKDB_CDN + "duckdb-browser-eh.worker.js",
      },
    });

    const worker = new Worker(bundle.mainWorker);
    const logger = new duckdb.ConsoleLogger();
    this._db = new duckdb.AsyncDuckDB(logger, worker);
    await this._db.instantiate(bundle.mainModule);

    this._conn = await this._db.connect();

    // Create a view over the remote Parquet file.
    // DuckDB handles HTTP range requests automatically.
    await this._conn.query(`
      CREATE VIEW annotations AS
      SELECT * FROM read_parquet('${this._parquetUrl}')
    `);
  }

  /**
   * Look up a specific variant by genomic coordinates.
   * @param {string} chrom - Chromosome (e.g. "12" or "chr12")
   * @param {number} pos   - 1-based position
   * @param {string} ref   - Reference allele
   * @param {string} alt   - Alternate allele
   * @returns {Promise<object[]>} Array of annotation rows
   */
  async lookupVariant(chrom, pos, ref, alt) {
    const normalizedChrom = chrom.replace(/^chr/i, "");
    const result = await this._conn.query(`
      SELECT * FROM annotations
      WHERE chrom = '${normalizedChrom}'
        AND pos = ${pos}
        AND ref = '${ref}'
        AND alt = '${alt}'
      ORDER BY is_canonical DESC, impact ASC
    `);
    return this._toObjects(result);
  }

  /**
   * Look up all variants for a gene.
   * @param {string} geneName - Gene symbol (e.g. "KRAS")
   * @returns {Promise<object[]>} Array of annotation rows
   */
  async lookupGene(geneName) {
    const result = await this._conn.query(`
      SELECT * FROM annotations
      WHERE gene_name = '${geneName}'
      ORDER BY chrom_numeric, pos
      LIMIT 1000
    `);
    return this._toObjects(result);
  }

  /**
   * Execute a free-form SQL query against the annotations view.
   * @param {string} sql - SQL query
   * @returns {Promise<object[]>} Array of result rows
   */
  async query(sql) {
    const result = await this._conn.query(sql);
    return this._toObjects(result);
  }

  /** Close the DuckDB connection and database. */
  async close() {
    if (this._conn) {
      await this._conn.close();
      this._conn = null;
    }
    if (this._db) {
      await this._db.terminate();
      this._db = null;
    }
  }

  /** Convert a DuckDB Arrow result to plain JS objects. */
  _toObjects(result) {
    const rows = [];
    for (let i = 0; i < result.numRows; i++) {
      const row = {};
      for (const field of result.schema.fields) {
        const col = result.getChild(field.name);
        row[field.name] = col.get(i);
      }
      rows.push(row);
    }
    return rows;
  }
}
