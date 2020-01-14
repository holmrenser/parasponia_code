#!/usr/bin/env gt

--[[
  Copyright (c) 2014-2015 Sascha Steinbiss <ss34@sanger.ac.uk>
  Copyright (c) 2014-2015 Genome Research Ltd

  Permission to use, copy, modify, and distribute this software for any
  purpose with or without fee is hereby granted, provided that the above
  copyright notice and this permission notice appear in all copies.

  THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES
  WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF
  MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR
  ANY SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
  WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN
  ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF
  OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.
]]

function usage()
  io.stderr:write("For a GFF3 file, output EMBL file suitable for Chado loader.\n")
  io.stderr:write(string.format("Usage: %s <GFF annotation> <organism name> <sequence>\n" , arg[0]))
  os.exit(1)
end

if #arg < 3 then
  usage()
end

function gff3_encode(s)
  s = string.gsub(s, "[\t\n\r;=%&,]", function (c)
        return string.format("%%%02X", string.byte(c))
      end)
  return s
end

function gff3_decode(s)
  if not s then
    return s
  end
  s = string.gsub(s, "%%([0-9a-fA-F][1-9a-fA-F])", function (n)
        return string.char(tonumber("0x" .. n))
      end)
  return s
end

function gff3_explode(tab)
  local ret = {}
  for _,item in ipairs(tab) do
    local _tmp = {}
    for k,v in pairs(item) do
      table.insert(_tmp, k .. "=" .. v)
    end
    table.insert(ret, gff3_encode(table.concat(_tmp, ";")))
  end
  return table.concat(ret, ",")
end

function gff3_extract_structure(str)
  ret = {}
  for _,v in ipairs(split(str, ", ?")) do
    res = {}
    v = gff3_decode(v)
    for _,pair in ipairs(split(v, ";")) do
      if string.len(pair) > 0 then
        key, value = unpack(split(pair, "="))
        res[key] = value
      end
    end
    table.insert(ret, res)
  end
  return ret
end

function file_exists(name)
   local f=io.open(name,"r")
   if f~=nil then io.close(f) return true else return false end
end

function split(str, pat)
   local t = {}
   local fpat = "(.-)" .. pat
   local last_end = 1
   local s, e, cap = str:find(fpat, 1)
   while s do
      if s ~= 1 or cap ~= "" then
        table.insert(t,cap)
      end
      last_end = e+1
      s, e, cap = str:find(fpat, last_end)
   end
   if last_end <= #str then
      cap = str:sub(last_end)
      table.insert(t, cap)
   end
   return t
end

function trim(s)
  -- from PiL2 20.4
  return (s:gsub("^%s*(.-)%s*$", "%1"))
end

function string:split(sep)
  local sep, fields = sep or ":", {}
  local pattern = string.format("([^%s]+)", sep)
  self:gsub(pattern, function(c) fields[#fields+1] = c end)
  return fields
end

function get_fasta(filename, sep)
  local keys = {}
  local seqs = {}
  local cur_hdr = nil
  local cur_seqs = {}
  if not file_exists(filename) then
    error("file " .. filename .. " can not be loaded")
  end
  for l in io.lines(filename) do
    hdr = l:match(">(.*)")
    if hdr then
      hdr = trim(split(hdr,"%s+")[1])
      table.insert(keys, hdr)
      if #cur_seqs > 0 and cur_hdr then
        if not seqs[cur_hdr] then
          seqs[cur_hdr] = table.concat(cur_seqs, sep)
        end
      end
      cur_hdr = hdr
      cur_seqs = {}
    else
      table.insert(cur_seqs, l)
    end
  end
  if cur_hdr and not seqs[cur_hdr] then
    seqs[cur_hdr] = table.concat(cur_seqs, sep)
  end
  return keys, seqs
end

function get_fasta_nosep(filename)
  return get_fasta(filename, "")
end

region_mapping = gt.region_mapping_new_seqfile_matchdescstart(arg[3])

collect_vis = gt.custom_visitor_new()
collect_vis.pps = {}
function collect_vis:visit_feature(fn)
  if fn:get_type() == "polypeptide" then
    local df = fn:get_attribute("Derives_from")
    if df then
      self.pps[df] = fn
    end
  end
  return 0
end

function format_embl_attrib(node, attrib, qualifier, fct)
  if node and node:get_attribute(attrib) then
    for _,v in ipairs(split(node:get_attribute(attrib), ",")) do
      if fct then
        s = fct(v)
      else
        s = gff3_decode(v)
      end
      if s then
        io.write("FT                   /" .. qualifier .. "=\"" .. s .."\"\n")
      end
    end
  end
end

function format_embl_sequence(sequence)
  local a = 0
  local c = 0
  local g = 0
  local t = 0
  local other = 0
  local l = 0
  -- count char distribution
  for ch in sequence:gmatch("%a") do
    ch = string.lower(ch)
    l = l + 1
    if ch == "a" then
      a = a + 1
    elseif ch == "c" then
      c = c + 1
    elseif ch == "g" then
      g = g + 1
    elseif ch == "t" then
      t = t + 1
    else
      other = other + 1
    end
  end
  -- show statistics
  io.write("SQ   Sequence " .. l .. " BP; " .. a .. " A; ".. c .. " C; "..
                                            g .. " G; ".. t .. " T; "..
                                            other .. " other;\n")
  local i = 1
  local pos = 0
  -- format and output sequence
  io.write("     ")
  for c in sequence:gmatch("%a", 10) do
    io.write(c)
    if i % 10 == 0 then
      io.write(" ")
    end
    if i % 60 == 0 then
      io.write(string.format("%9s\n     ", i))
    end
    i = i + 1
  end
  io.write(string.format(string.rep(' ',(80-i%60-(i%60)/10-13)) .. "%10d\n", i-1))
end

embl_vis = gt.custom_visitor_new()
embl_vis.pps = collect_vis.pps
embl_vis.last_seqid = nil
function embl_vis:visit_feature(fn)
  if embl_vis.last_seqid ~= fn:get_seqid() then
    if embl_vis.last_seqid then
      format_embl_sequence(collect_vis.seqs[embl_vis.last_seqid])
      io.write("//\n")
      io.output():close()
    end
    embl_vis.last_seqid = fn:get_seqid()
    io.output(fn:get_seqid()..".embl", "w+")
    io.write("ID   " .. fn:get_seqid() .. "; XXX; linear; XXX; XXX; XXX; XXX.\n")
    io.write("XX   \n")
    io.write("DE   " .. arg[2] .. ", " .. fn:get_seqid() .. "\n")
    io.write("XX   \n")
    io.write("AC   XXX.\n")
    io.write("XX   \n")
    io.write("PR   Project:00000000;\n")
    io.write("XX   \n")
    io.write("KW   \n")
    io.write("XX   \n")
    io.write("RN   [1]\n")
    io.write("RA   Authors;\n")
    io.write("RT   Title;\n")
    io.write("RL   Unpublished.\n")
    io.write("XX   \n")
    io.write("OS   " .. arg[2] .. "\n")
    io.write("XX   \n")
    io.write("FH   Key             Location/Qualifiers\n")
    io.write("FH   \n")
    io.write("FT   source          1.." .. collect_vis.lengths[fn:get_seqid()] .. "\n")
    io.write("FT                   /organism=\"" .. arg[2] .. "\"\n")
    io.write("FT                   /mol_type=\"genomic DNA\"\n")
  end

  for node in fn:get_children() do
    if node:get_type() == "mRNA" or node:get_type() == "pseudogenic_transcript" then
      local cnt = 0
      for cds in node:get_children() do
        if cds:get_type() == "CDS" or cds:get_type() == "pseudogenic_exon" then
          cnt = cnt + 1
        end
      end
      io.write("FT   CDS             ")
      if node:get_strand() == "-" then
        io.write("complement(")
      end
      if cnt > 1 then
        io.write("join(")
      end
      local i = 1
      local coding_length = 0
      for cds in node:get_children() do
        if cds:get_type() == "CDS" or cds:get_type() == "pseudogenic_exon" then
          if i == 1 and fn:get_attribute("Start_range") then
            io.write("<")
          end
          io.write(cds:get_range():get_start())
          io.write("..")
          if i == cnt and fn:get_attribute("End_range") then
            io.write(">")
          end
          io.write(cds:get_range():get_end())
          if i ~= cnt then
            io.write(",")
          end
          coding_length = coding_length + cds:get_range():length()
          i = i + 1
        end
      end
      if cnt > 1 then
        io.write(")")
      end
      if node:get_strand() == "-" then
        io.write(")")
      end
      io.write("\n")
      local pp = self.pps[node:get_attribute("ID")]
      format_embl_attrib(pp, "product", "product",
          function (s)
            local pr_a = gff3_extract_structure(s)
            local gprod = pr_a[1].term
            if gprod then
              return gprod
            else
              return nil
            end
          end)
      format_embl_attrib(node , "ID", "locus_tag", nil)
      if fn:get_type() == "pseudogene" then
        io.write("FT                   /pseudo\n")
        --io.write("FT                   /pseudogene=\"unknown\"\n")
      end
      format_embl_attrib(pp, "Dbxref", "EC_number",
          function (s)
            m = string.match(s, "EC:([0-9.-]+)")
            if m then
              return m
            else
              return nil
            end
          end)
      -- add gene to 'unroll' multiple spliceforms
      local geneid = fn:get_attribute("ID")
      if geneid then
        io.write("FT                   /gene=\"".. geneid .. "\"\n")
      end
      -- translation
      local protseq = nil
      if node:get_type() == "mRNA" then
        protseq = node:extract_and_translate_sequence("CDS", true,
                                                      region_mapping)
        io.write("FT                   /translation=\"" .. protseq:sub(1,-2) .."\"\n")
      end
      io.write("FT                   /transl_table=1\n")
      -- add name
      local name = fn:get_attribute("Name")
      if name then
        io.write("FT                   /gene_synonym=\"".. name .. "\"\n")
      end
    elseif node:get_type() == "tRNA" then
      io.write("FT   tRNA            ")
      if node:get_strand() == "-" then
        io.write("complement(")
      end
      io.write(node:get_range():get_start() .. ".." .. node:get_range():get_end())
      if node:get_strand() == "-" then
        io.write(")")
      end
      io.write("\n")
      if node:get_attribute("aa") then
        io.write("FT                   /product=\"" .. node:get_attribute("aa") .. " transfer RNA")
        if node:get_attribute("anticodon") then
          io.write(" (" .. node:get_attribute("anticodon") .. ")")
        end
        io.write("\"\n")
      end
      format_embl_attrib(node , "ID", "locus_tag", nil)
    elseif string.match(node:get_type(), "snRNA") or string.match(node:get_type(), "snoRNA") then
      io.write("FT   ncRNA            ")
      if node:get_strand() == "-" then
        io.write("complement(")
      end
      io.write(node:get_range():get_start() .. ".." .. node:get_range():get_end())
      if node:get_strand() == "-" then
        io.write(")")
      end
      io.write("\n")
      io.write("FT                   /ncRNA_class=\"" .. node:get_type() .. "\"\n")
      format_embl_attrib(node , "ID", "locus_tag", nil)
    elseif string.match(node:get_type(), "rRNA") then
      io.write("FT   rRNA            ")
      if node:get_strand() == "-" then
        io.write("complement(")
      end
      io.write(node:get_range():get_start() .. ".." .. node:get_range():get_end())
      if node:get_strand() == "-" then
        io.write(")")
      end
      io.write("\n")
      io.write("FT                   /product=\"" .. node:get_type() .. "\"\n")
      format_embl_attrib(node , "ID", "locus_tag", nil)
    end
  end
  return 0
end

-- make simple visitor stream, just applies given visitor to every node
visitor_stream = gt.custom_stream_new_unsorted()
function visitor_stream:next_tree()
  local node = self.instream:next_tree()
  if node then
    node:accept(self.vis)
  end
  return node
end

-- get sequences, sequence lengths and derives_from relationships
visitor_stream.instream = gt.gff3_in_stream_new_sorted(arg[1])
visitor_stream.vis = collect_vis
local gn = visitor_stream:next_tree()
while (gn) do
  gn = visitor_stream:next_tree()
end
keys, seqs = get_fasta_nosep(arg[3])
collect_vis.seqs = {}
collect_vis.lengths = {}
for k,v in pairs(seqs) do
  collect_vis.seqs[k] = v
  collect_vis.lengths[k] = v:len()
end

-- output EMBL code as we traverse the GFF
visitor_stream.instream = gt.gff3_in_stream_new_sorted(arg[1])
visitor_stream.vis = embl_vis
embl_vis.gaf = gaf_store
local gn = visitor_stream:next_tree()
while (gn) do
  gn = visitor_stream:next_tree()
end
-- output last seq
format_embl_sequence(collect_vis.seqs[embl_vis.last_seqid])
io.write("//\n")
io.output():close()
