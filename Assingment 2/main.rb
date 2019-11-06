require './assingment2_objects/search_int_net'
require 'rest-client'
require 'json'
require 'ruby-progressbar'

# check number of input files and if name are duplicated
if ARGV.uniq.length != 2
  abort "\nError: wrong number of arguments.\n\nPlease use the command: ruby main.rb gene_codes_list.txt report.txt"
end

puts "\nWelcome to Assingment 2!"

puts "\nChecking files..."
agi_file = ARGV[0]
puts agi_file
fnameout = ARGV[1]
puts 'Succes!'

puts "\nThis program will find protein-protein interactions between genes from
your list and will predict possible interaction networks.It will return a report
file with all the interaction networks found and member will be annotated with
KEGG and GO functional annotations.\n"

puts "\nYou can choose the IntAct MIscore cutoff. By default is 0.485.
Do you want to change it?\n"

ans = $stdin.gets.chomp
if Regexp.new(/[Yy]es|[Ss]i|[Ss]|[Yy]/).match(ans)
  puts "\nWhich cutoff would you like to stablished? (try 0.37 to see a few examples of direct interactions).\n"
  cutoff = $stdin.gets.chomp.to_f
  puts 'This may be longer but more interactions will be found' if cutoff < 0.485
  puts 'This may be faster but less interactions will be found' if cutoff > 0.485
else
  cutoff = 0.485
end
puts "\nStarting program..."
interactors = Assingment2Objects::InteractionNetwork.search_interactors(agi_file)
interactions = Assingment2Objects::InteractionNetwork.search_ppi(interactors, cutoff)
pos_net = Assingment2Objects::InteractionNetwork.possible_networks(interactions)

# START REPORT
File.open(fnameout, 'w+') do |f|
  f.puts '--------------------------------------------------------------------------------------'
  f.puts "  Total interaction networks found: #{pos_net.length} (deepest level = 3) (PPI scores >= #{cutoff})"
  f.puts '--------------------------------------------------------------------------------------'
  f.puts '-> Direct networks are those in which all members belong to the gene list.'
  f.puts "-> Indirect networks are those in which members A and C belong to the gene list\n   but not B."
  f.puts ''
  f.puts "Direct Networks: #{pos_net.find_all { |net| '(direct)' == net.predicted_net.match(/\(\w+\)/).to_s }.length}"
  f.puts "Indirect Networks: #{pos_net.find_all { |net| '(indirect)' == net.predicted_net.match(/\(\w+\)/).to_s }.length}"
  f.puts ''
  (pos_net.each { |net| f.puts "Network #{net.network_num}: #{net.predicted_net}" })
  pos_net.each do |net|
    puts ''
    f.puts '----------------------------------------------------------------------------------'
    f.puts "                   Network #{net.network_num} (deep level = #{net.members.length})"
    f.puts '----------------------------------------------------------------------------------'
    f.puts net.predicted_net.to_s
    f.puts '-----------'
    f.puts '| Members |'
    f.puts '-----------'
    f.puts net.members.map(&:agi_code).join(' ').to_s
    f.puts '---------------------------'
    f.puts '| Annotations for network |'
    f.puts '---------------------------'
    f.puts "\nGO terms in common with all members:"
    if net.go_common.empty?
      f.puts 'Common GO annotations not found.'
    else
      f.puts("GO_ID\tGo_Term")
      net.go_common.each { |id, term| f.puts "#{id.strip}\t#{term.strip}" }
    end
    f.puts "\nTotal KEGG terms:"
    if net.kegg_all.empty?
      f.puts 'KEGG annotations not found.'
    else
      f.puts("KEGG_ID\tPathway_Name")
      net.kegg_all.each { |id, term| f.puts "#{id.strip}\t#{term.strip}" }
    end
    f.puts ''
    f.puts '---------------------------'
    f.puts '| Annotations for members |'
    f.puts '---------------------------'
    net.members.each do |member|
      f.puts "Gene ID: #{member.agi_code}"
      if member.kegg.nil? || member.kegg.empty?
        f.puts('KEGG annotation not found')
      else
        f.puts("KEGG_ID\tPathway_Name")
        member.kegg.each { |id, name| f.puts "#{id.strip}\t#{name.strip}" }
      end
      if member.go.nil? || member.go.empty?
        f.puts('GO annotation not found')
      else
        f.puts("GO_ID\tGo_Term")
        member.go.each { |id, term| f.puts "#{id.strip}\t#{term.strip}" }
      end
      f.puts ''
    end
  end
end
# ENDS REPORT

puts "\nReport saved in #{Dir.pwd}\n"
puts "\nEnd of program."
